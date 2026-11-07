"""A searchable, non-destructive auxiliary-software log and energy summary."""
import argparse
from datetime import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import re
import tempfile

NAME='总流程.log'
SKIP_DIRS={'scratch','OPTIM','__pycache__','.git'}
QUANTUM_STAGES={'opt','opt1','opt2','opt3','freq','sp','coupling','coupling_shielding'}
TEXT_SUFFIXES={'.log','.out','.csv','.txt','.json','.dat','.mdp','.chg','.top','.itp','.tcl'}
COORD_SUFFIXES={'.xyz','.xtc','.trr','.gro','.pdb','.mol2','.sdf','.edr','.tpr','.cpt'}


def quantum_file(path):
    """Recognize legacy quantum names too; helper logs remain auxiliary output."""
    path=Path(path);name=path.name.lower()
    if name.endswith(('_multiwfn.log','_formchk.log')):return False
    if name.startswith(('gau_','orca_')):return True
    stage=path.parent.name.lower()
    if stage in QUANTUM_STAGES or stage.startswith(('nmr','ecd')):return True
    if stage in ('resp_gas','resp_solvent') and name in ('gau.out','job.out','job.log'):return True
    return path.suffix.lower() in {'.gjf','.chk','.fchk','.gbw','.rwf','.hess','.wfn','.wfx'}


def _atomic(path,text):
    path.parent.mkdir(parents=True,exist_ok=True)
    fd,name=tempfile.mkstemp(prefix='.'+path.name+'.',dir=str(path.parent))
    try:
        with os.fdopen(fd,'w',encoding='utf-8',newline='\n') as f:f.write(text)
        os.replace(name,path)
    finally:
        if os.path.exists(name):os.unlink(name)


def _json(path):
    try:return json.loads(path.read_text(encoding='utf-8-sig'))
    except (OSError,ValueError):return None


def ensemble_table(path,limit):
    """Parse full ensembles strictly; never call a truncated frame complete."""
    energies=[]
    with path.open(encoding='utf-8') as f:
        while True:
            header=f.readline()
            if not header:break
            if not header.strip():continue
            n=int(header);comment=f.readline()
            if n<1 or not comment:raise ValueError('不完整XYZ头')
            match=re.match(r'\s*(?:Energy\s*=\s*)?([-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?)',comment)
            if match is None:raise ValueError('注释中没有可识别的Hartree能量')
            energy=float(match[1].replace('D','E').replace('d','e'))
            if not math.isfinite(energy):raise ValueError('非有限能量')
            for _ in range(n):
                row=f.readline().split()
                if len(row)!=4:raise ValueError('XYZ帧正在写入或不完整')
                if not all(math.isfinite(float(v)) for v in row[1:]):raise ValueError('非有限坐标')
            energies.append(energy)
    if not energies:raise ValueError('空集合')
    ranked=sorted(enumerate(energies,1),key=lambda item:(item[1],item[0]))
    low=ranked[0][1]
    lines=['构象数={}；电子能排序，不含热校正；相对能单位 kcal/mol。'.format(len(energies)),
           '排名\t原集合ID\tE/Hartree\t相对E/kcal mol-1\t本预算内']
    for rank,(index,energy) in enumerate(ranked,1):
        lines.append('{}\tc{:04d}\t{:.12f}\t{:.6f}\t{}'.format(rank,index,energy,(energy-low)*627.509474,
                     '是' if limit==0 or rank<=limit else '否'))
    return '\n'.join(lines)


def build(project,cfg=None,out=None,live=None):
    project=Path(project).resolve()
    if cfg is None:
        config=next((project/name for name in ('job.json','job_compute.json') if (project/name).is_file()),None)
        if config is None:raise ValueError('找不到项目配置')
        cfg=json.loads(config.read_text(encoding='utf-8'))
    out=Path(out if out is not None else project/cfg['output']).resolve()
    target=project/'运行日志'/NAME
    resources=project/'workflow/resources'
    candidates=cfg.get('candidates',[])
    files={};inventory=[];tables=[];warnings=[];quantum=[]
    def add(path,group,static=False):
        if not path.is_file() or path.is_symlink():return
        if any(part in SKIP_DIRS or part.startswith('.') for part in path.relative_to(out).parts):return
        if quantum_file(path):
            if path.suffix.lower() in ('.log','.out'):quantum.append(path)
            return
        suffix=path.suffix.lower()
        if suffix in COORD_SUFFIXES:
            inventory.append((path,group))
            return
        if suffix not in TEXT_SUFFIXES:return
        if static and (resources/path.name).is_file() and path.suffix!='.mdp':return
        files[str(path)]=(path,group)
    for candidate in candidates:
        cid=candidate['id'];run=out/cid
        for folder,label in [('sampling','前处理、动力学、GFN搜索与聚类'),('dft_dedup','DFT后聚类'),('shermo','热校正与自由能分布')]:
            root=run/folder
            if root.is_dir():
                for p in sorted(root.rglob('*')):add(p,cid+' / '+label,static=folder=='sampling')
        for folder in sorted((run/'conformers').glob('c*')):
            if not folder.is_dir():continue
            for stage in folder.iterdir():
                if not stage.is_dir():continue
                for p in stage.glob('*'):
                    if p.is_file() and quantum_file(p) and p.suffix.lower() in ('.log','.out'):quantum.append(p)
        for name in ('state.json','conformers.json'):
            p=run/name
            if p.is_file():inventory.append((p,cid+' / 计算记录'))
        for stage in ('gfn1','gfn2'):
            p=run/'sampling'/stage/'cluster.xyz'
            if p.is_file():
                try:tables.append((cid+' / '+stage.upper(),p,ensemble_table(p,cfg['max_conformers'])))
                except (OSError,ValueError,OverflowError) as exc:warnings.append(str(p)+'：'+str(exc))
        manifest=_json(run/'conformers.json')
        if manifest:
            rows=manifest.get('records',[])
            text=['仅针对已完成频率验收的所算集合；不是全部构象空间覆盖率。',
                  '构象\t优化电子能/Hartree\tORCA电子能/Hartree\tG/Hartree\t布居/%']
            for row in rows:text.append('\t'.join(str(v) for v in (row['id'],row.get('opt_energy',''),row.get('sp_energy',''),row.get('gibbs_hartree',''),100*row['weight'])))
            tables.append((cid+' / 最终自由能比例',run/'conformers.json','\n'.join(text)))
    for root in (out/'exports',out/'module_exports/conformers'):
        if root.is_dir():
            for path in sorted(root.rglob('*')):add(path,'能量与构象分布导出')
    events=project/'运行日志/software_events.jsonl'
    lines=['Auto ConfoNMR 总流程日志',
           '更新时间：'+datetime.now().isoformat(timespec='seconds'),
           '项目：'+str(project),
           '范围：动力学、CREST、聚类、能量表、Shermo与布居等辅助过程的原始文本。',
           'Gaussian/ORCA正文保留在独立文件；本日志仅列出其位置。轨迹与结构数据保留原文件。',
           '这是当前可用输出的汇总快照，正在写入的阶段可能尚未完成。',
           '\n[计算条件]',json.dumps({k:cfg.get(k) for k in ('protocol','charge','multiplicity','max_conformers','clustering_profile','clustering','md_seeds','shermo','modules')},ensure_ascii=False,indent=2)]
    if live:
        lines+=['\n[当前阶段]']
        for row in live.get('candidates',[]):
            lines.append('{}：{}；当前={}；DFT优化完成={}；频率合格={}；ORCA完成={}'.format(
                row['id'],row['status'],row['current'],row['optimized'],row['positive_frequencies'],row['single_points']))
    lines+=['\n[软件调用时间与返回码]']
    if events.is_file():
        try:lines.append(events.read_text(encoding='utf-8',errors='replace'))
        except OSError as exc:warnings.append(str(exc))
    else:lines.append('尚无软件调用事件。')
    lines+=['\n[原始输出正文索引]']
    for number,(path,group) in enumerate(files.values(),1):lines.append('{}\t{}\t{}'.format(number,group,path))
    for number,(path,group) in enumerate(files.values(),1):
        try:
            # Read a byte snapshot; preserve all content without arbitrary tail limits.
            raw=path.read_bytes()
            body=raw.decode('utf-8-sig',errors='replace')
            lines+=['\n'+'='*88,'[原始输出 {}] {}'.format(number,group),'文件：'+str(path),
                    '快照字节数={}；SHA256={}'.format(len(raw),hashlib.sha256(raw).hexdigest()),'-'*88,body,
                    '[原始输出 {} 结束]'.format(number)]
        except OSError as exc:warnings.append(str(path)+'：'+str(exc))
    for label,path,text in tables:lines+=['\n'+'='*88,'[能量/分布表] '+label,'来源：'+str(path),text]
    lines+=['\n[Gaussian/ORCA独立原始输出索引：不合并正文]']
    lines += [str(p) for p in sorted(set(quantum))]
    lines+=['\n[轨迹、结构及计算记录：保留原文件]']
    for path,group in inventory:
        try:lines.append('{}\t{}\t{} bytes'.format(group,path,path.stat().st_size))
        except OSError:pass
    if warnings:lines+=['\n[暂不可读取/未完成]',*warnings]
    text='\n'.join(lines)+'\n'
    _atomic(target,text)
    return {'path':str(target),'auxiliary_text_files':len(files),'quantum_log_links':len(set(quantum)),
            'energy_tables':len(tables),'bytes':len(text.encode('utf-8')),'warnings':warnings}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('project')
    args=parser.parse_args();print(json.dumps(build(args.project),ensure_ascii=False))
