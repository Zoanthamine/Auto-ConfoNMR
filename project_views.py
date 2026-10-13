#!/usr/bin/env python3
"""Publish readable, non-destructive views of a project's original calculation files."""
import argparse
import csv
import html
import json
import os
from pathlib import Path
import re
import shutil
import tempfile
import time
from urllib.parse import quote
from progress_status import snapshot


def atomic_text(path, text):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    fd,name=tempfile.mkstemp(prefix='.'+path.name+'.',dir=str(path.parent))
    try:
        with os.fdopen(fd,'w',encoding='utf-8') as f:f.write(text)
        os.replace(name,path)
    finally:
        if os.path.exists(name):os.unlink(name)


def load_config(project):
    project=Path(project).resolve()
    p=next((project/name for name in ['job.json','job_compute.json'] if (project/name).is_file()),None)
    if p is None:raise ValueError('找不到 job.json 或 job_compute.json：'+str(project))
    c=json.loads(p.read_text(encoding='utf-8'))
    return p,c,(project/c['output']).resolve()


def link_file(source,target,copy=False):
    """Do not modify or remove any original calculation file."""
    source=Path(source).resolve();target=Path(target)
    if not source.is_file():return '不存在'
    target.parent.mkdir(parents=True,exist_ok=True)
    if copy:
        if target.exists() and not target.is_symlink():
            a,b=source.stat(),target.stat()
            if a.st_size==b.st_size and a.st_mtime_ns==b.st_mtime_ns:return '复制快照'
        fd,name=tempfile.mkstemp(prefix='.copy_',dir=str(target.parent));os.close(fd)
        try:
            shutil.copy2(source,name);os.replace(name,target)
        finally:
            if os.path.exists(name):os.unlink(name)
    elif target.is_symlink():
        if target.resolve()!=source:raise ValueError('查看链接指向另一批计算：'+str(target))
    elif target.exists():
        # An earlier finalized copy is retained, never replaced with a live link.
        return '保留快照（运行04刷新）'
    else:
        target.symlink_to(source)
    return '复制快照' if copy else '原文件实时链接'


def geometry(log):
    txt=Path(log).read_text(errors='replace')
    if 'Normal termination of Gaussian' not in txt or 'Error termination' in txt:return None
    blocks=re.findall(r'(?:Standard|Input) orientation:\s*\n\s*-+\n.*?\n\s*-+\n(.*?)\n\s*-+',txt,re.S)
    if not blocks:return None
    symbols='X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr'.split()
    atoms=[]
    for line in blocks[-1].splitlines():
        p=line.split()
        if len(p)!=6 or int(p[1])>=len(symbols):return None
        atoms.append(symbols[int(p[1])]+' '+' '.join(p[3:6]))
    return str(len(atoms))+'\nExtracted from '+str(log)+'\n'+'\n'.join(atoms)+'\n'


def event(project,kind,**values):
    try:
        p=Path(project)/'运行日志'/'software_events.jsonl';p.parent.mkdir(parents=True,exist_ok=True)
        with p.open('a',encoding='utf-8') as f:
            f.write(json.dumps(dict(time=time.strftime('%Y-%m-%d %H:%M:%S'),event=kind,**values),ensure_ascii=False)+'\n')
    except OSError as exc:
        print('软件事件索引写入失败，独立原始日志仍保留：'+str(exc),flush=True)


def collect(project,copy_logs=False):
    project=Path(project).resolve();_,cfg,out=load_config(project)
    live=snapshot(project,cfg,out)
    live_by_id={r['id']:r for r in live['candidates']}
    root=project/'过程文件';root.mkdir(exist_ok=True)
    direct = cfg.get('output_layout') == 'direct_v1' and out == root.resolve()
    index=[];weights=[];summary=[]
    from protocol_options import method_map
    stages={'opt':'01_优化','opt1':'01_优化_尝试1','opt2':'01_优化_尝试2','opt3':'01_优化_尝试3','freq':'02_频率'}
    stages.update({'nmr'+mid:'{:02d}_NMR_{}_{}'.format(int(mid)+2,mid,f) for mid,(f,_,_) in method_map(cfg).items()})
    light={'.log','.out','.gjf','.inp','.json','.txt','.csv','.mdp','.tcl','.top','.itp','.chg','.xyz','.pdb','.gro'}
    def expose(src,dst,label):
        if direct:
            src=Path(src).resolve()
            if root.resolve() in src.parents and root.resolve() in Path(dst).resolve().parents:
                index.append([str(src),str(src.relative_to(project)),label,'实际计算文件',src.stat().st_size])
            else:
                form=link_file(src,dst,True)
                index.append([str(src),str(dst.relative_to(project)),label,form,src.stat().st_size])
            return
        # Trajectories/ensemble XYZ can be large; leave live links even when finalizing.
        copy=copy_logs and src.suffix.lower() in light and src.stat().st_size<=32*1024*1024
        form=link_file(src,dst,copy)
        index.append([str(src),str(dst.relative_to(project)),label,form,src.stat().st_size])
    for c in cfg['candidates']:
        cid=c['id']
        if re.fullmatch(r'[A-Za-z0-9_-]+',cid) is None:raise ValueError('候选名称不安全')
        run=out/cid;sampling=run/'sampling'
        source=(project/c['input']).resolve()
        if source.is_file():expose(source,root/'03_构型'/cid/('Input_'+cid+source.suffix),'输入结构')
        for phase in ['resp_gas','resp_solvent']:
            d=sampling/phase
            if d.exists():
                for f in d.iterdir():
                    if f.is_file() and f.suffix in ('.gjf','.out','.log','.chk','.fchk','.json','.chg','.in'):
                        expose(f,root/'01_Gaussian'/cid/'00_RESP'/phase/f.name,'RESP '+phase)
        for d in sorted((run/'conformers').glob('c*')):
            if not d.is_dir():continue
            for stage,label in stages.items():
                sd=d/stage
                if not sd.is_dir():continue
                for f in sd.iterdir():
                    if f.is_file() and f.name!='.workflow.lock':expose(f,root/'01_Gaussian'/cid/d.name/label/f.name,label)
                logs = list(sd.glob('Gau_*.log')) if direct else [sd/'job.log']
                if copy_logs and logs and logs[0].is_file():
                    xyz=geometry(logs[0])
                    if xyz:
                        target=root/'03_构型'/cid/('Geom_{}_{}_{}.xyz'.format(stage,cid,d.name.upper()));atomic_text(target,xyz)
                        index.append([str(logs[0]),str(target.relative_to(project)),label+'末次几何','从正常终止日志提取XYZ',len(xyz.encode())])
            sd=d/'sp'
            if sd.is_dir():
                for f in sd.iterdir():
                    if f.is_file():expose(f,root/'02_ORCA'/cid/d.name/f.name,'ORCA单点')
        if sampling.is_dir():
            for f in sampling.rglob('*'):
                if f.is_file() and (f.suffix in light or f.name in ('receipt.json',)):
                    if f.name in ('trajectory.xyz','trajectory_stereo_checked.xyz','crest_ensemble.xyz','stereo_checked.xyz','cluster.xyz','gfn2.xyz'):
                        expose(f,root/'03_构型'/cid/'采样集合'/f.relative_to(sampling),'采样集合')
                    elif f.suffix in ('.log','.json','.mdp','.top','.itp','.chg','.tcl'):
                        expose(f,root/'04_采样与软件日志'/cid/f.relative_to(sampling),'采样/软件记录')
        for f in (run/'shermo').glob('*'):
            if f.is_file():expose(f,root/'05_Shermo与分布'/cid/f.name,'Shermo')
        for name in ['conformers.json','state.json']:
            f=run/name
            if f.is_file():expose(f,root/'05_Shermo与分布'/cid/name,name)
        m=run/'conformers.json';s=run/'state.json'
        state=json.loads(s.read_text())['status'] if s.is_file() else '尚未开始'
        if m.is_file():
            z=json.loads(m.read_text());records=z['records'];g0=min(r['gibbs_hartree'] for r in records)
            summary.append([cid,state,z['ensemble_total'],z['dft_selected'],len(records),len(z['nmr_ids']),100*z['coverage_within_computed_ensemble']])
            for r in records:weights.append([cid,r['id'],(r['gibbs_hartree']-g0)*627.509474,100*r['weight'],100*r['nmr_weight'] if 'nmr_weight' in r else '',r['id'] in z['nmr_ids']])
        else:
            r=live_by_id[cid]
            summary.append([cid,state,r['gfn2'] if r['gfn2'] is not None else '待生成',
                            r['dft_inputs'],r['positive_frequencies'],r['nmr_complete_conformers'],'待最终汇总'])
    for f in (out/'exports').glob('*'):
        if f.is_file():expose(f,project/'最终结果'/'最新导出'/f.name,'当前exports')
    # Preserve the case's previously published upper/lower comparisons as named views.
    for point in ['upper','lower']:
        for f in (project/'results'/point).glob('*'):
            if f.is_file():expose(f,project/'最终结果'/('上分离点' if point=='upper' else '下分离点')/f.name,point)
    if copy_logs:
        from experiment_batch import collect_comparisons
        comparisons=collect_comparisons(project)
        for row in comparisons:
            for f in (project/row['result']).rglob('*'):
                if f.is_file():
                    index.append([str(f),str(f.relative_to(project)),'实验比较 '+row['sample'],'复制快照',f.stat().st_size])
    if direct:
        # Index software-generated fixed names too, without copies or symlinks.
        known={row[0] for row in index}
        for f in sorted(root.rglob('*')):
            if f.is_file() and not f.is_symlink() and str(f.resolve()) not in known and 'scratch' not in f.parts:
                index.append([str(f),str(f.relative_to(project)),'过程记录','实际计算文件',f.stat().st_size])
    def csv_text(headers,rows):
        import io
        s=io.StringIO(newline='');w=csv.writer(s);w.writerow(headers);w.writerows(rows)
        return '\ufeff'+s.getvalue()
    atom_rows=[]
    for candidate in cfg['candidates']:
        mapping=out/candidate['id']/'sampling/atom_map.json'
        if mapping.is_file():
            for row in json.loads(mapping.read_text(encoding='utf-8')):
                atom_rows.append([candidate['id'],row['calculation_atom'],row['element'],row['input_serial']])
    atomic_text(project/'原子编号对照.csv',csv_text(['candidate','calculation_atom','element','input_serial'],atom_rows))
    atomic_text(project/'文件索引.csv',csv_text(['原始文件','集中查看位置','阶段','形式','字节数'],index))
    atomic_text(project/'构象分布与混合比例.csv',csv_text(['结构','构象','相对G_kcal_mol','原始分布_percent','NMR混合_percent','用于NMR'],weights))
    storage_note = '本目录直接保存实际计算文件，无需依赖外部scratch链接。\n' if direct else '历史布局：查看目录不改变原始文件；大文件保持链接，普通输入/日志在整理时复制快照。\n'
    atomic_text(root/'原始计算位置.txt','原始计算目录：'+str(out)+'\n'+storage_note+'整理时间：'+time.strftime('%Y-%m-%d %H:%M:%S')+'\n')
    def ht(headers,rows):
        return '<table><thead><tr>'+''.join('<th>'+html.escape(str(x))+'</th>' for x in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+html.escape(str(x))+'</td>' for x in r)+'</tr>' for r in rows)+'</tbody></table>'
    links=''.join('<li><a href="'+quote(row[1])+'">'+html.escape(row[1])+'</a></li>' for row in index)
    page='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>NMR计算总览</title><style>body{font:16px "Microsoft YaHei",sans-serif;max-width:1200px;margin:32px auto;color:#243746;padding:0 24px}table{border-collapse:collapse;width:100%;margin:18px 0}th,td{border:1px solid #cbd6de;padding:9px;text-align:left}th{background:#174a67;color:white}a{color:#176a97}li{margin:6px 0}p{line-height:1.7}details{margin:24px 0}</style><h1>NMR 计算总览</h1>'''
    page+='<p>项目：'+html.escape(str(project))+'<br>更新：'+time.strftime('%Y-%m-%d %H:%M:%S')+'。状态来自已有记录；是否仍在运行请同时查询队列。</p>'
    page+=ht(['结构','状态','GFN2构象','DFT输入','正频率','NMR构象','覆盖/%'],summary)
    page+='<h2>阶段实时记录</h2><p>读取时间：'+html.escape(live['checked_at'])+'。只统计已完成回执；完整覆盖率等待最终汇总，运行状态同时以队列为准。</p>'
    page+=ht(['结构','当前阶段/构象','优化完成','ORCA完成','NMR01/02/03已完成任务','最近日志更新'],
             [[r['id'],r['current'],r['optimized'],r['single_points'],
               ' / '.join(str(r['nmr'].get(m,0)) for m in ('01','02','03')),r['log_modified'] or '尚无记录']
              for r in live['candidates']])
    page+='<p><a href="实时进度.html">02入口最近一次读取的详细进度快照</a></p>'
    page+='<h2>构象混合比例</h2><p>每候选内部混合；01/02/03及各实验样品共享相同分布。99%只针对已计算的DFT集合。</p>'+ht(['结构','构象','ΔG kcal/mol','原始/%','NMR/%','选中'],[[*r[:2],f'{r[2]:.5f}',f'{r[3]:.6f}',f'{r[4]:.6f}' if r[4]!='' else '未做NMR',r[5]] for r in weights])
    page+='<p><a href="文件索引.csv">完整文件索引</a> · <a href="构象分布与混合比例.csv">权重CSV</a></p><details><summary>展开各软件输入、输出和构型文件</summary><ul>'+links+'</ul></details>'
    page+='<p>力场警告请检查04_采样与软件日志中的topology_quality.json；正常结束与有限预算覆盖率不能证明采样充分。</p></html>'
    atomic_text(project/'计算总览.html',page)
    return {'files':len(index),'candidates':summary,'project':str(project)}


def refresh_from_env():
    project=os.environ.get('NMR_PROJECT_DIR')
    if not project:return
    try:collect(project,False)
    except Exception as exc:print('查看目录更新未完成（原始计算日志仍保留）：'+str(exc),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('project');p.add_argument('--copy',action='store_true')
    a=p.parse_args();print(json.dumps(collect(a.project,a.copy),ensure_ascii=False))
