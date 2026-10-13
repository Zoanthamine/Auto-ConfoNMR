#!/usr/bin/env python3
"""Chinese project settings and short manual entry points; science stays in nmrflow."""
import argparse
import configparser
from datetime import datetime
import getpass
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import shlex
from cluster_config import load_site, install_resources, runtime_text, missing_resources

from project_views import atomic_text, collect, load_config
from protocol_options import apply_choices, overview, LAYOUT, NMR_PROFILES, CLUSTERS

ROOT=Path(__file__).resolve().parent
PY=sys.executable
HOME=Path.home()
ENTRIES={'01_submit.sh':'submit','02_status.sh':'status','03_compare.sh':'compare','04_collect.sh':'collect'}


def settings_text(cfg):
    return '''# NMR中文配置；用文本编辑器修改并保存为UTF-8。
# 候选列表固定在job.json；新增/删除候选请建立新项目。
# 开始计算后，电荷、预算、输出目录等计算条件不能直接修改后续跑。
[方案]
# 粗分类coarse：GFN1 0.5/0.5；GFN2及DFT 0.3/0.2；T300。
# 细分类fine=GFN1 0.2/0.2；GFN2及DFT 0.2/0.1；T350。单位与全表见计算方案.md。
# 两档都不自动增加DFT预算；预算小可能让细分类的覆盖收益受限。
clustering_profile = {clustering_profile}
# legacy_pcm_g09=IEFPCM+G09Defaults（默认）；
# pcm_g16=PCM+G16默认值；smd_g16=SMD+G16默认值。
# 这里只改变NMR屏蔽步骤；不改变气相优化/频率/ORCA及Shermo。
nmr_profile = {nmr_profile}

[计算]
# 分子总电荷，中性为0；自旋多重度，闭壳层通常为1。
charge = {charge}
multiplicity = {multiplicity}
# 最多送去DFT的低GFN2能量构象数；0=全部，成本可能明显上升。
max_conformers = {max_conformers}
# 对已完成DFT集合的累计分布覆盖；0.99不代表整个构象空间覆盖99%。
population_coverage = {population_coverage}
# 新项目直接写入本项目“过程文件”，不是链接；无需修改此项。
# Gaussian/ORCA输入、完整输出、检查点、采样及权重都保存在这个真实目录中。
output = {output}
# 核数同步用于PBS单节点申请；Gaussian为总内存，ORCA为每进程内存。
# 下面的数值要符合节点和队列限制，ORCA maxcore×cores也需留有系统余量。
cores = {cores}
md_threads = {md_threads}
memory_gb = {memory_gb}
orca_maxcore_mb = {orca_maxcore_mb}
# Shermo分布温度/K和压力/atm。低频熵插值阈值100cm-1及因子1保留job.json默认。
temperature = {T}
pressure = {P}

[对比]
# 01=B3LYP/6-31G(d)→DP4+；02=mPW1PW91/6-311+G(d,p)→DP4+；03=wB97XD/6-31G(d)→STS。
# true=沿用当前自定义方案输出探索比较；false=导出输入并提示校准条件待审核。
# 所选方案的实际几何/分布/溶剂/版本差异会写入报告，不是标准校准的绝对可信度。
allow_exploratory = true
# 每个候选必须为每个实验signal提供原子归属；同一甲基3H放一个信号。
assignments = assignments.csv
# 自动处理根目录全部experiment_样品名.csv，不限制样品数量。
# 优先使用assignments_样品名.csv，找不到时才使用共用assignments.csv。
# 删除/移走不参与比较的实验表；不要把同一数据重复保存为不同文件。
experiment_pattern = experiment_*.csv
# 兼容单表experiment.csv；不存在时不读取。一般多样品使用上面的命名即可。
experiment = experiment.csv
'''.format(**cfg,T=cfg['shermo']['T'],P=cfg['shermo']['P'])


def install_entries(project,workflow,python=None):
    for name,action in ENTRIES.items():
        script='''#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\\n\\t'
PROJECT_DIR="$(cd -- "$(dirname -- "${{BASH_SOURCE[0]}}")" && pwd)"
readonly PROJECT_DIR
exec {python} "$PROJECT_DIR/workflow/project_console.py" {action} "$PROJECT_DIR" "$@"
'''.format(python=shlex.quote(python or PY),action=action)
        note = {'compare':'# 自动检查所有实验表后逐一比较；仅检查用：bash 03_compare.sh --check-only\n', 'collect':'# 收集所有已完成样品及历史比较结果，不提交计算。\n'}.get(action, '')
        script=script.replace('set -Eeuo pipefail',note+'set -Eeuo pipefail')
        atomic_text(project/name,script);(project/name).chmod(0o755)


def manifest(project,workflow):
    files=[p for p in workflow.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='submissions.tsv']
    atomic_text(project/'程序版本SHA256.json',json.dumps({str(p.relative_to(workflow)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},ensure_ascii=False,indent=2)+'\n')


def prepare(project,structures,charge=0,multiplicity=1,cap=6,nmr_profile='legacy_pcm_g09',clustering_profile='coarse'):
    if nmr_profile not in NMR_PROFILES or clustering_profile not in CLUSTERS:raise ValueError('无效方案，尚未建立项目')
    project=Path(project).resolve()
    if project==ROOT.resolve() or ROOT.resolve() in project.parents:
        raise ValueError('项目必须建在程序目录外，防止递归复制；建议建在个人主目录下')
    if project.exists():raise ValueError('项目目录已经存在，请换一个新名字：'+str(project))
    structures=[Path(p).resolve() for p in structures]
    if not structures or not all(p.is_file() for p in structures):raise ValueError('未找到明确指定的输入结构')
    site=load_site()
    missing=missing_resources(site)
    if missing:raise ValueError('外部资源缺少：'+', '.join(missing))
    subprocess.run([sys.executable,str(ROOT/'init_project.py'),str(project),*[str(p) for p in structures],
                    '--charge',str(charge),'--multiplicity',str(multiplicity)],check=True)
    wf=project/'workflow'
    shutil.copytree(ROOT,wf,ignore=shutil.ignore_patterns('__pycache__','*.pyc','legacy_original','submissions.tsv','.git','.site-resources','.history'))
    install_resources(site,wf)
    atomic_text(wf/'config/cluster.json',json.dumps(site,ensure_ascii=False,indent=2)+'\n')
    atomic_text(wf/'cluster_runtime.sh',runtime_text(site))
    p=wf/'run_nmr.pbs';s=p.read_text()
    pattern=r"(?m)^readonly WORKFLOW_HOME='[^'\n]*'$"
    if len(re.findall(pattern,s))!=1:raise ValueError('PBS程序路径格式已改变，请保留项目并检查')
    atomic_text(p,re.sub(pattern,lambda _:"readonly WORKFLOW_HOME='"+str(wf)+"'",s))
    cfg=json.loads((project/'job.json').read_text())
    cfg['output']='过程文件';cfg['output_layout']='direct_v1';cfg['max_conformers']=cap
    cfg['software']=dict(site['software'])
    apply_choices(cfg,nmr_profile,clustering_profile)
    for b in cfg['postprocessing']:b['allow_exploratory']=True
    atomic_text(project/'job.json',json.dumps(cfg,ensure_ascii=False,indent=2)+'\n')
    atomic_text(project/'project.ini',settings_text(cfg))
    atomic_text(project/'计算方案.md',overview(cfg))
    preview_inputs(project,cfg,wf)
    install_entries(project,wf,site['python']);manifest(project,wf)
    guide=ROOT/'README.md'
    if guide.exists():shutil.copy2(guide,project/'使用说明.md')
    collect(project,False)
    print('\nAuto ConfoNMR 项目已经准备好：'+str(project)+'\n请编辑 project.ini，上传实验CSV；然后运行 01_submit.sh。没有提交计算。')


def sync(project):
    path,cfg,out=load_config(project)
    if path.name!='job.json':raise ValueError('这是已有完成项目的计算配置；仅支持查看/整理，不能用新入口覆盖后重新提交')
    ini=Path(project)/'project.ini'
    if not ini.is_file():ini=Path(project)/'项目配置.ini'  # Read-only naming compatibility for tests and earlier snapshots.
    if not ini.is_file():raise ValueError('缺少项目配置；请先通过00_new_project.sh建立项目')
    p=configparser.ConfigParser(interpolation=None);p.read(ini,encoding='utf-8-sig')
    known_calc={'charge','multiplicity','max_conformers','population_coverage','output','cores','md_threads','memory_gb','orca_maxcore_mb','temperature','pressure'}
    known_post={'allow_exploratory','assignments','experiment','experiment_pattern','upper_experiment','lower_experiment','upper_assignments','lower_assignments'}
    required_post={'allow_exploratory','assignments','experiment'}
    if set(p.sections())!={'方案','计算','对比'} or set(p['方案'])!={'clustering_profile','nmr_profile'} or set(p['计算'])!=known_calc or not required_post <= set(p['对比']) or not set(p['对比']) <= known_post:
        raise ValueError('配置的节名或字段有遗漏/拼写错误，请参考原始project.ini')
    if cfg.get('method_layout')!=LAYOUT:raise ValueError('该项目使用历史编号；请用原固定程序，或建立新项目')
    apply_choices(cfg,p.get('方案','nmr_profile').strip(),p.get('方案','clustering_profile').strip())
    for key in ['charge','multiplicity','max_conformers','cores','md_threads','memory_gb','orca_maxcore_mb']:cfg[key]=p.getint('计算',key)
    cfg['population_coverage']=p.getfloat('计算','population_coverage');cfg['output']=p.get('计算','output').strip()
    if not cfg['output']:raise ValueError('output不能为空')
    if cfg.get('output_layout')=='direct_v1' and (Path(project)/cfg['output']).resolve() != (Path(project)/'过程文件').resolve():
        raise ValueError('本版原始输出固定为本项目的过程文件目录；更换存储位置请调整新项目根目录')
    cfg['shermo']['T']=p.getfloat('计算','temperature');cfg['shermo']['P']=p.getfloat('计算','pressure')
    cfg['experiment']=p.get('对比','experiment').strip();cfg['assignments']=p.get('对比','assignments').strip()
    for b in cfg['postprocessing']:b['allow_exploratory']=p.getboolean('对比','allow_exploratory')
    from nmrflow import validate_config,calculation_config
    validate_config(cfg)
    resolved=json.loads(json.dumps(cfg));resolved['output']=str((Path(project)/cfg['output']).resolve())
    for c in resolved['candidates']:c['input']=str((Path(project)/c['input']).resolve())
    # Do not lose a previous identity merely by editing output to a different directory.
    prior=json.loads(path.read_text());oldout=(Path(project)/prior['output']).resolve()
    for c in prior['candidates']:
        state=oldout/c['id']/'state.json'
        if state.exists():
            existing=json.loads(state.read_text())['identity']['config']
            if existing!=calculation_config(resolved):
                raise ValueError('已有计算记录，计算配置发生变化。请新建项目；只改实验表/对比开关不需要新计算。')
    atomic_text(path,json.dumps(cfg,ensure_ascii=False,indent=2)+'\n')
    atomic_text(Path(project)/'计算方案.md',overview(cfg))
    preview_inputs(Path(project),cfg,ROOT)
    return path,cfg,p


def preview_inputs(project,cfg,workflow):
    """Publish exact NMR route choices, without inventing molecular coordinates."""
    from protocol_options import nmr_text
    for mid in cfg['methods']:
        source=Path(workflow)/'resources'/('06_tem'+mid+'_nmr.gjf')
        text=nmr_text(source.read_text(),cfg)
        text=re.sub(r'%mem=.*','%mem={}GB'.format(cfg['memory_gb']),text,flags=re.I)
        text=re.sub(r'%nprocshared=.*','%nprocshared={}'.format(cfg['cores']),text,flags=re.I)
        text=re.sub(r'(?m)^0 1$','{} {}'.format(cfg['charge'],cfg['multiplicity']),text)
        atomic_text(Path(project)/'参数预览'/('NMR_'+mid+'_模板_含坐标占位符.gjf'),text)


def queue_guard(project):
    """Fail closed on scheduler outages and reject an active recorded submission."""
    result=subprocess.run(['qstat','-u',getpass.getuser()],capture_output=True,text=True)
    if result.returncode:raise ValueError('调度服务无法核实，暂不提交，防止重复作业：'+result.stderr.strip())
    record=Path(project)/'提交记录.json'
    if record.exists():
        for row in json.loads(record.read_text()):
            q=subprocess.run(['qstat','-f',row['job_id']],capture_output=True,text=True)
            if q.returncode==0:
                m=re.search(r'job_state\s*=\s*(\w+)',q.stdout)
                if not m or m[1] not in ('C','F'):raise ValueError('已有作业活动中，拒绝重复提交：'+row['job_id'])
            elif not re.search(r'Unknown Job|Unknown job|Invalid Job',q.stderr+q.stdout,re.I):
                raise ValueError('无法确认旧作业状态，暂不提交：'+row['job_id'])


def submit(project):
    path,cfg,_=sync(project)
    _,_,out=load_config(project)
    complete=all((out/c['id']/'state.json').is_file() and json.loads((out/c['id']/'state.json').read_text())['status']=='complete' for c in cfg['candidates'])
    if complete:raise ValueError('本项目已经计算完成。复看用02，改实验数据用03，不重新提交。')
    queue_guard(project)
    print('检查配置、软件和输入，通过后提交一次……',flush=True)
    r=subprocess.run(['bash',str(ROOT/'submit.sh'),str(path)],cwd=str(project),capture_output=True,text=True)
    stamp=datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    atomic_text(Path(project)/'运行日志'/('submit_'+stamp+'.log'),r.stdout+r.stderr)
    print(r.stdout,end='');print(r.stderr,end='',file=sys.stderr)
    if r.returncode:raise ValueError('检查或提交未成功；先查看运行日志/submit文件，不盲目再次提交')
    m=re.search(r'Submitted\s+(\S+)\.',r.stdout)
    if not m:raise ValueError('提交返回值无法识别；请检查submissions.tsv和队列，禁止直接重提')
    record=Path(project)/'提交记录.json';rows=json.loads(record.read_text()) if record.exists() else []
    rows.append({'job_id':m[1],'submitted':stamp,'config_sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    atomic_text(record,json.dumps(rows,ensure_ascii=False,indent=2)+'\n')
    print('以后运行02_status.sh查看；进度见运行日志/run_*.log；Gaussian/ORCA完整输出见过程文件中的独立日志。')


def status(project):
    collect(project,False)
    from progress_status import main as live_status
    live_status([str(project),'--html'])
    print('打开计算总览.html，或进入过程文件查看；这不会提交计算。')


def compare(project,check_only=False):
    _,cfg,p=sync(project)
    if not check_only:
        queue_guard(project)
    from experiment_batch import run_batch
    result=run_batch(project,cfg,p,check_only=check_only)
    if not check_only:
        collect(project,True)
        print('本批实验比较已保存；请核对各样品的数据与结果复核.md及信号误差.csv。')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare','submit','status','compare','collect']);p.add_argument('project',nargs='?');p.add_argument('--inputs',nargs='+');p.add_argument('--charge',type=int,default=0);p.add_argument('--multiplicity',type=int,default=1);p.add_argument('--cap',type=int,default=6)
    p.add_argument('--nmr-profile',choices=list(NMR_PROFILES),default='legacy_pcm_g09')
    p.add_argument('--clustering-profile',choices=list(CLUSTERS),default='coarse')
    p.add_argument('--check-only',action='store_true',help='仅检查实验表与原子归属（用于compare）')
    a=p.parse_args()
    if a.action=='prepare':
        if a.project and a.inputs:return prepare(a.project,a.inputs,a.charge,a.multiplicity,a.cap,a.nmr_profile,a.clustering_profile)
        project=input('新项目英文名称（例如NMR_Molecule_20261013）：').strip()
        if not re.fullmatch(r'[A-Za-z0-9_-]+',project):raise ValueError('项目名称只能含英文、数字、下划线或短横线')
        folder=Path(input('已上传PDB所在目录：').strip()).expanduser();files=sorted(folder.glob('*.pdb'))
        if not files:raise ValueError('目录中没有.pdb文件')
        print('将只使用以下结构：\n'+'\n'.join(str(x) for x in files))
        if input('确认这些是全部候选？输入yes：').strip().lower()!='yes':raise ValueError('未确认，不建立项目')
        charge=int(input('总电荷[默认0]：') or 0);multi=int(input('自旋多重度[默认1]：') or 1);cap=int(input('DFT最多构象数[默认6，0全部]：') or 6)
        choice=input('聚类档位：1粗分类（默认）；2细分类：').strip() or '1'
        if choice not in ('1','2'):raise ValueError('聚类选项应为1或2')
        nmr=input('NMR方案：1 IEFPCM+G09兼容（默认）；2 PCM+G16；3 SMD+G16探索：').strip() or '1'
        if nmr not in ('1','2','3'):raise ValueError('NMR选项应为1、2或3')
        prepare(Path(load_site()['projects_root'])/project,files,charge,multi,cap,{'1':'legacy_pcm_g09','2':'pcm_g16','3':'smd_g16'}[nmr],{'1':'coarse','2':'fine'}[choice])
    else:
        if not a.project:p.error('需要项目目录')
        project=Path(a.project).resolve()
        if a.action=='submit':submit(project)
        elif a.action=='status':status(project)
        elif a.action=='compare':compare(project,a.check_only)
        else:print(json.dumps(collect(project,True),ensure_ascii=False))


if __name__=='__main__':
    try:main()
    except (ValueError,OSError,subprocess.SubprocessError) as exc:sys.exit('停止：'+str(exc))
