"""Discover experimental samples, validate explicit assignments and retain independent runs."""
import contextlib
import csv
from datetime import datetime
import io
import json
import math
from pathlib import Path
import re
import shutil

from project_views import atomic_text
from nmrflow import digest, output_lock

EXPERIMENT_COLUMNS = ['signal', 'nucleus', 'ppm']
ASSIGNMENT_COLUMNS = ['candidate', 'signal', 'atoms', 'hybridization', 'carbon_type', 'heavy_atom']


def read_table(path, columns):
    with Path(path).open(encoding='utf-8-sig', newline='') as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None or len(reader.fieldnames) != len(columns) or set(reader.fieldnames) != set(columns):
            raise ValueError(str(path)+': 表头必须为 '+','.join(columns))
        rows = []
        for line, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError('{}:{}: 列数与表头不符'.format(path, line))
            if any(value != value.strip() for value in row.values()):
                raise ValueError('{}:{}: 单元格含首尾空格，请删除后重新保存'.format(path,line))
            rows.append(row)
    if not rows:
        raise ValueError(str(path)+': 数据表为空，请填写真实数据，或移走不用的实验表')
    return rows


def discover(project, settings, allow_empty=False):
    project = Path(project)
    pattern = settings.get('对比', 'experiment_pattern', fallback='experiment_*.csv').strip()
    if not pattern or '/' in pattern or '\\' in pattern or '..' in pattern:
        raise ValueError('experiment_pattern 必须为项目根目录的文件名模式，如 experiment_*.csv')
    found = {}
    # Module CSVs are separate schemas, never NMR shift samples (even when disabled).
    reserved={(project/name).resolve() for name in ('experiment_j.csv','experiment_ecd.csv')}
    for key in ('j_experiment','ecd_experiment'):
        name=settings.get('模块后处理',key,fallback='').strip()
        if name:reserved.add((project/name).resolve())
    def add(path, sample, assignment=''):
        path = path.resolve()
        if not path.is_file():
            return
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*', sample):
            raise ValueError('实验样品名只使用英文、数字、下划线或短横线：'+sample)
        for old in found.values():
            if old['experiment'] == path:
                return
        if sample.casefold() in found:
            raise ValueError('实验样品名冲突（含大小写）：'+sample)
        specific = project/('assignments_'+sample+'.csv')
        mapping = project/assignment if assignment else specific if specific.is_file() else project/settings.get('对比', 'assignments')
        if not mapping.is_file():
            raise ValueError(sample+': 缺少原子归属表 '+str(mapping))
        found[sample.casefold()] = {'sample': sample, 'experiment': path, 'assignments': mapping.resolve()}
    # Explicit legacy names are supported, but no longer limit the sample count.
    for key, sample, akey in [('experiment', 'single', 'assignments'),
                             ('upper_experiment', 'upper', 'upper_assignments'),
                             ('lower_experiment', 'lower', 'lower_assignments')]:
        name = settings.get('对比', key, fallback='').strip()
        if name:
            add(project/name, sample, settings.get('对比', akey, fallback='').strip() if key != 'experiment' else '')
    for path in sorted(project.glob(pattern)):
        if path.resolve() in reserved or path.name.endswith(".template.csv"):continue
        if path.name.startswith('experiment_') and path.suffix == '.csv':
            add(path, path.stem[len('experiment_'):])
        else:
            raise ValueError('实验文件应命名为 experiment_样品名.csv：'+path.name)
    if not found and not allow_empty:
        raise ValueError('没有实验数据。填写 experiment_样品名.csv 和对应 assignments_样品名.csv（或共用 assignments.csv）。')
    return sorted(found.values(), key=lambda item: item['sample'].casefold())


def atom_elements(out, candidate):
    mapping = Path(out)/candidate/'sampling/atom_map.json'
    if mapping.is_file():
        records = json.loads(mapping.read_text(encoding='utf-8'))
        result = {r['calculation_atom']: r['element'] for r in records}
        if len(result) != len(records) or set(result) != set(range(1, len(result)+1)):
            raise ValueError(candidate+': atom_map 原子编号不连续或重复')
        return result
    manifest = Path(out)/candidate/'conformers.json'
    if manifest.is_file():
        geometry = json.loads(manifest.read_text(encoding='utf-8'))['records'][0]['geometry']
        return {i: atom[0] for i, atom in enumerate(geometry, 1)}
    raise ValueError(candidate+': 尚无原子映射/构象记录，请等待结构准备完成后再检查数据')


def validate_sample(item, cfg, out):
    from postprocess import validate_branch_data
    exp = read_table(item['experiment'], EXPERIMENT_COLUMNS)
    mapping = read_table(item['assignments'], ASSIGNMENT_COLUMNS)
    signals = {row['signal']: row for row in exp}
    if len(signals) != len(exp) or '' in signals:
        raise ValueError('signal 不能为空或重复')
    warnings = []
    for row in exp:
        if row['nucleus'] not in ('C', 'H'):
            raise ValueError(row['signal']+': nucleus 只填写 C 或 H')
        try:
            ppm = float(row['ppm'])
        except ValueError as exc:
            raise ValueError(row['signal']+': ppm 必须为数字，不填写单位或区间') from exc
        if not math.isfinite(ppm):
            raise ValueError(row['signal']+': ppm 不能为 NaN/Inf')
        lo, hi = (-20, 250) if row['nucleus'] == 'C' else (-5, 20)
        if not lo <= ppm <= hi:
            warnings.append(row['signal']+': 位移超出常见范围，请核对核种、单位、小数点和实验来源')
    for nucleus in ('C', 'H'):
        values = [float(row['ppm']) for row in exp if row['nucleus'] == nucleus]
        if values and (len(values) < 3 or max(values) == min(values)):
            raise ValueError(nucleus+': 每种参与拟合的核至少三条信号，且位移不能全部相同')
    candidates = [c['id'] for c in cfg['candidates']]
    keys = [(r['candidate'], r['signal']) for r in mapping]
    expected = {(candidate, signal) for candidate in candidates for signal in signals}
    if len(keys) != len(set(keys)):
        raise ValueError('同一候选与 signal 的归属行重复')
    if set(keys) != expected:
        raise ValueError('归属表应恰好覆盖每个候选的全部实验信号；缺少={}；多余={}'.format(
            sorted(expected-set(keys)), sorted(set(keys)-expected)))
    assembled, checked = {}, []
    for candidate in candidates:
        elements = atom_elements(out, candidate)
        used = set()
        assembled[candidate] = []
        for row in (r for r in mapping if r['candidate'] == candidate):
            if not re.fullmatch(r'[1-9][0-9]*(;[1-9][0-9]*)*', row['atoms']):
                raise ValueError(candidate+':'+row['signal']+': atoms 使用从1开始的整数，等价原子用分号分隔')
            atoms = [int(a) for a in row['atoms'].split(';')]
            if len(set(atoms)) != len(atoms) or used & set(atoms):
                raise ValueError(candidate+':'+row['signal']+': 原子重复或被分配到多个信号')
            used.update(atoms)
            nucleus = signals[row['signal']]['nucleus']
            if any(elements.get(atom) != nucleus for atom in atoms):
                raise ValueError(candidate+':'+row['signal']+': 原子序号不存在或核种与元素不一致')
            if len(atoms) > 1:
                warnings.append(candidate+':'+row['signal']+': 等价原子组 '+row['atoms']+' 需要人工确认化学等价性')
            assembled[candidate].append(dict(row, nucleus=nucleus))
            checked.append(dict(row, nucleus=nucleus, ppm=signals[row['signal']]['ppm']))
        omitted = sorted(a for a, el in elements.items() if el in ('C', 'H') and a not in used)
        if omitted:
            warnings.append(candidate+': 未纳入比较的 C/H 原子 '+','.join(map(str, omitted))+'；确认缺失信号、重叠峰或有意排除的原因')
    parameters = json.loads((Path(__file__).parent/'calibration.json').read_text())
    for branch in cfg['postprocessing']:
        validate_branch_data(branch, assembled, parameters)
        if branch['kind'] in ('dp4', 'dp4plus') and len(candidates) < 2:
            raise ValueError('DP4+ 相对比较至少需要两个候选；一个实验样品也可比较多个候选')
        if branch['kind'] == 'sts' and (len(candidates) > 12 or sum(e['nucleus'] == 'C' for e in exp) > 91):
            raise ValueError('上游 STS 工作簿最多12个候选、91条碳信号；实验样品文件数量没有此限制')
        if branch['kind'] == 'dp4plus' and (len(candidates) > 128 or len(exp) > 300):
            raise ValueError('上游 DP4+ 工作簿容量超限')
    for signal in signals:
        rows = [r for r in mapping if r['signal'] == signal]
        fields = ['hybridization'] + (['carbon_type', 'heavy_atom'] if signals[signal]['nucleus'] == 'C' else [])
        if any(len({r[field] for r in rows}) != 1 for field in fields):
            raise ValueError(signal+': 上游工作簿共用分类，候选间分类不一致；请审核真实分类，不要为通过检查改分类')
    return {'sample': item['sample'], 'experiment': str(item['experiment']), 'assignments': str(item['assignments']),
            'experiment_sha256': digest(item['experiment']), 'assignments_sha256': digest(item['assignments']),
            'atom_map_sha256': {c: digest(Path(out)/c/'sampling/atom_map.json') for c in candidates if (Path(out)/c/'sampling/atom_map.json').is_file()},
            'signals': len(exp), 'candidates': candidates, 'warnings': warnings, 'checked_assignments': checked}


def csv_text(columns, rows):
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=columns)
    writer.writeheader(); writer.writerows(rows)
    return '\ufeff'+stream.getvalue()


def review_results(destination, report):
    """Expose fitted signal residuals, while keeping all statistical formulas unchanged."""
    from postprocess import table
    result = json.loads((destination/'comparison_results.json').read_text())
    residuals = []
    warnings = list(report['warnings'])
    for method, candidates in result['descriptive_fits'].items():
        inputs = table(destination/('method'+method+'_comparison_input.csv'))
        for candidate, nuclei in candidates.items():
            for nucleus, fit in nuclei.items():
                if fit['slope'] >= 0:
                    warnings.append('method'+method+'/'+candidate+'/'+nucleus+': 屏蔽-位移拟合斜率非负，请核对归属和输入是否误填位移')
                for row in (r for r in inputs if r['nucleus'] == nucleus):
                    fitted = fit['slope']*float(row[candidate])+fit['intercept']
                    residuals.append({'method': method, 'candidate': candidate, 'signal': row['signal'], 'nucleus': nucleus,
                                      'experimental_ppm': float(row['experimental_ppm']), 'fitted_ppm': fitted,
                                      'residual_ppm': fitted-float(row['experimental_ppm'])})
    atomic_text(destination/'信号误差.csv', csv_text(
        ['method', 'candidate', 'signal', 'nucleus', 'experimental_ppm', 'fitted_ppm', 'residual_ppm'], residuals))
    lines = ['# 实验数据与比较结果复核', '', '样品：'+report['sample'], '',
             '此清单和信号误差用于发现问题，不代替实验归属审核。信号误差是对当前数据拟合后的残差，不能视作外部预测误差。', '',
             '## 需要人工核对', '',
             '- 核对实验样品名称、溶剂、温度、ppm 单位及原始谱图，确认未混入其他产物的信号。',
             '- 对照 atom_map.json 和结构，逐项核对计算原子号、非等价氢、等价原子组与 DP4+/STS 分类。',
             '- 查看 method01/02/03_comparison_input.csv，确认每个候选使用相同信号集合。',
             '- 对照信号误差.csv 查找异常峰；任何排除或重新归属均需实验依据并保留理由。',
             '- 核对两组 DP4+ 与 STS 的概率、误差及条件差异；排序冲突时先查输入与计算质量。',
             '- 概率是在给定候选集合内的相对评价，不代表样品纯度或异构体浓度；高概率不能排除未计算的结构。',
             '- Excel 打开后重算公式，并与 comparison_results.json 核对；不要仅依据工作簿缓存值。', '',
             '## 自动检查提示', '']
    lines += ['- '+warning for warning in warnings] or ['- 未发现上述格式问题；化学归属仍需人工核对。']
    lines += ['', '## 拟合残差较大的信号（仅用于排查）', '', '|方法|候选|信号|核种|残差/ppm|', '|---|---|---|---|---|']
    for row in sorted(residuals, key=lambda r: abs(r['residual_ppm']), reverse=True)[:20]:
        lines.append('|{method}|{candidate}|{signal}|{nucleus}|{residual_ppm:.5f}|'.format(**row))
    atomic_text(destination/'数据与结果复核.md', '\n'.join(lines)+'\n')


def run_batch(project, cfg, settings, check_only=False):
    from postprocess import export
    project = Path(project).resolve()
    resolved = json.loads(json.dumps(cfg))
    out = (project/cfg['output']).resolve()
    resolved['output'] = str(out)
    for candidate in resolved['candidates']:
        candidate['input'] = str((project/candidate['input']).resolve())
    items = discover(project, settings)
    reports, errors = [], []
    for item in items:
        try:
            reports.append(validate_sample(item, resolved, out))
        except (ValueError, KeyError, OSError) as exc:
            errors.append(item['sample']+': '+str(exc))
    stamp = datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    check = {'checked_at': stamp, 'status': 'invalid' if errors else 'format_and_mapping_checked',
             'samples': reports, 'errors': errors, 'note': '自动检查不证明化学归属正确；比较后还需人工复核。'}
    atomic_text(project/'运行日志'/('data_check_'+stamp+'.json'), json.dumps(check, ensure_ascii=False, indent=2)+'\n')
    if errors:
        raise ValueError('实验表检查未通过，未进行本批比较：\n'+'\n'.join(errors))
    for item, report in zip(items, reports):
        print('样品 {}：{}；归属 {}；{} 条信号；{} 项人工复核提示'.format(
            item['sample'], item['experiment'].name, item['assignments'].name, report['signals'], len(report['warnings'])), flush=True)
    if check_only:
        print('检查完成，未执行统计比较。详见运行日志/data_check_'+stamp+'.json')
        return check
    with output_lock(out):
        for candidate in cfg['candidates']:
            state = out/candidate['id']/'state.json'
            if not state.is_file() or json.loads(state.read_text())['status'] != 'complete':
                raise ValueError(candidate['id']+': 计算未完成，拒绝进行本批比较')
        run = project/'results'/('comparison_'+stamp)
        run.mkdir(parents=True, exist_ok=False)
        index = {'run': run.name, 'status': 'running', 'samples': []}
        index_path = run/'batch.json'
        def record():
            atomic_text(index_path, json.dumps(index, ensure_ascii=False, indent=2)+'\n')
        record()
        try:
            for item, report in zip(items, reports):
                name = item['sample']; folder = run/name; inputs = run/('_inputs_'+name)
                inputs.mkdir()
                for key in ('experiment', 'assignments'):
                    shutil.copy2(item[key], inputs/(key+'.csv'))
                    if digest(inputs/(key+'.csv')) != report[key+'_sha256']:
                        raise ValueError(name+': 检查后输入发生变化，请重新运行03')
                point = json.loads(json.dumps(resolved))
                point.update(experiment=str(inputs/'experiment.csv'), assignments=str(inputs/'assignments.csv'))
                atomic_text(inputs/'job_post.json', json.dumps(point, ensure_ascii=False, indent=2)+'\n')
                atomic_text(inputs/'data_check.json', json.dumps(report, ensure_ascii=False, indent=2)+'\n')
                log = project/'运行日志'/('compare_'+name+'_'+stamp+'.log')
                entry = {'sample': name, 'status': 'running', 'result': str(folder.relative_to(project)),
                         'inputs': str(inputs.relative_to(project)), 'log': str(log.relative_to(project))}
                index['samples'].append(entry); record()
                try:
                    with log.open('w', encoding='utf-8') as stream, contextlib.redirect_stdout(stream), contextlib.redirect_stderr(stream):
                        export(point, out, project, destination=folder)
                        review_results(folder, report)
                    entry['status'] = 'complete'; record()
                    atomic_text(project/'最终结果'/('最新_'+name+'_结果位置.txt'), str(folder)+'\n')
                    print(name+' 完成：'+str(folder)+'；请阅读 数据与结果复核.md', flush=True)
                except Exception as exc:
                    entry.update(status='failed', error=str(exc)); record()
                    with log.open('a',encoding='utf-8') as stream:
                        stream.write('\n比较失败：'+str(exc)+'\n')
                    raise
            index['status'] = 'complete'; record()
            atomic_text(project/'results/latest_batch.json', json.dumps(index, ensure_ascii=False, indent=2)+'\n')
        except Exception as exc:
            index.update(status='failed', error=str(exc)); record()
            raise
    return index


def collect_comparisons(project):
    """Copy every completed sample in every run; preserve history and index latest per sample."""
    project = Path(project).resolve()
    rows = []
    latest = {}
    for manifest in sorted((project/'results').glob('comparison_*/batch.json')):
        index = json.loads(manifest.read_text(encoding='utf-8'))
        for entry in index['samples']:
            if entry['status'] != 'complete':
                continue
            source = (project/entry['result']).resolve()
            if project/'results' not in source.parents or not source.is_dir():
                raise ValueError('比较结果路径无效：'+str(source))
            target = project/'最终结果/实验比较'/index['run']/entry['sample']
            if source.is_symlink():
                raise ValueError('比较结果不能为链接：'+str(source))
            for file in source.iterdir():
                if file.is_file():
                    from project_views import link_file
                    link_file(file, target/file.name, copy=True)
            inputs = project/entry['inputs']
            for file in inputs.iterdir():
                if file.is_file():
                    from project_views import link_file
                    link_file(file, target/'输入记录'/file.name, copy=True)
            rows.append({'sample': entry['sample'], 'run': index['run'], 'batch_status': index['status'],
                         'result': str(target.relative_to(project)), 'original': entry['result']})
            latest[entry['sample']] = rows[-1]
    atomic_text(project/'最终结果/实验比较索引.csv', csv_text(['sample', 'run', 'batch_status', 'result', 'original'], rows))
    atomic_text(project/'最终结果/各样品最新结果.csv', csv_text(['sample', 'run', 'batch_status', 'result', 'original'], list(latest.values())))
    return rows
