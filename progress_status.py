#!/usr/bin/env python3
"""Read-only live progress. Completion receipts are evidence, never a new calculation."""
import argparse
from datetime import datetime
import json
import html
import os
from pathlib import Path
import re
import subprocess
import tempfile


def read_json(path,warnings):
    if not path.is_file():return None
    try:return json.loads(path.read_text(encoding='utf-8'))
    except (OSError,ValueError) as exc:
        warnings.append('记录暂不可读：{}（{}）'.format(path.name,type(exc).__name__))
        return None


def xyz_count(path,warnings):
    """A growing/truncated XYZ is unknown, not a silently shortened ensemble."""
    try:
        count=0
        with path.open(encoding='utf-8') as f:
            while True:
                line=f.readline()
                if not line:break
                if not line.strip():continue
                n=int(line)
                if n<1 or not f.readline():raise ValueError('XYZ header')
                for _ in range(n):
                    if len(f.readline().split())!=4:raise ValueError('XYZ atom')
                count+=1
        return count
    except (OSError,ValueError):
        warnings.append('构象文件正在写入或暂不可读：'+str(path))
        return None


def recent_events(project):
    f=project/'运行日志/software_events.jsonl'
    if not f.is_file():return []
    try:
        with f.open('rb') as h:
            h.seek(max(0,f.stat().st_size-131072))
            lines=h.read().decode('utf-8',errors='replace').splitlines()
    except OSError:return []
    rows=[]
    for line in lines:
        try:row=json.loads(line)
        except ValueError:continue
        if isinstance(row,dict):rows.append(row)
    return rows


def stage_label(cwd,run):
    parts=cwd.relative_to(run).parts
    stage=parts[-1] if parts else ''
    labels={'freq':'频率','sp':'ORCA单点','resp_gas':'气相RESP','resp_solvent':'溶剂RESP',
            'gfn1':'GFN1优化/聚类','gfn2':'GFN2优化/聚类','coupling':'J耦合','coupling_shielding':'J分支屏蔽'}
    if re.fullmatch(r'nmr\d+',stage):label=stage.upper()
    elif re.fullmatch(r'opt\d*',stage):label='DFT优化'+stage[3:]
    elif re.fullmatch(r'ecd\d+',stage):label=stage.upper()
    else:label=labels.get(stage,stage or '阶段切换')
    conformer=next((v for v in parts if re.fullmatch(r'c\d+',v)),None)
    return (conformer+' / ' if conformer else '')+label


def candidate_snapshot(candidate,cfg,output,events):
    run=output/candidate['id'];warnings=[]
    state=read_json(run/'state.json',warnings)
    status=state.get('status','状态记录暂不可读') if isinstance(state,dict) else ('状态记录暂不可读' if (run/'state.json').exists() else '尚未开始')
    manifest=read_json(run/'conformers.json',warnings)
    if not isinstance(manifest,dict):manifest=None
    stages={};inputs=set();rejected=set()
    for d in sorted((run/'conformers').glob('c*')):
        if not d.is_dir() or not re.fullmatch(r'c\d+',d.name):continue
        if any(d.glob('source.xyz')) or any(d.glob('Geom_source_*.xyz')):inputs.add(d.name)
        for stage in d.iterdir():
            if not stage.is_dir():continue
            if stage.name.startswith('opt') and (any(stage.glob('*.gjf')) or (stage/'receipt.json').is_file()):inputs.add(d.name)
            receipt=read_json(stage/'receipt.json',warnings)
            if not isinstance(receipt,dict):continue
            result=receipt.get('status','accepted')
            if result=='rejected_frequency':rejected.add(d.name)
            if result not in ('accepted','complete','success'):continue
            key='opt' if re.fullmatch(r'opt\d*',stage.name) else stage.name
            stages.setdefault(key,set()).add(d.name)
    gfn2=manifest.get('ensemble_total') if manifest else None
    if gfn2 is None:
        ensemble=next((f for f in (run/'sampling/gfn2.xyz',run/'sampling/gfn2/cluster.xyz') if f.is_file()),None)
        if ensemble:gfn2=xyz_count(ensemble,warnings)
    methods=cfg.get('methods',['01','02','03'])
    nmr={m:len(stages.get('nmr'+m,set())) for m in methods}
    completed=set.intersection(*(stages.get('nmr'+m,set()) for m in methods)) if methods else set()
    target=len(manifest.get('nmr_ids',[]))*len(methods) if manifest else None
    coverage=manifest.get('coverage_within_computed_ensemble') if manifest else None
    current='全部计算完成' if status=='complete' else ('尚未开始' if status=='尚未开始' else '阶段切换/汇总中（请结合队列）')
    event=None
    for row in reversed(events):
        cwd=Path(row.get('cwd',''))
        try:cwd.relative_to(run)
        except ValueError:continue
        event=row;break
    log=None;modified=None
    if event:
        if status!='complete':
            current=stage_label(Path(event['cwd']),run)
            current+='（运行中）' if event.get('event')=='started' else '（已结束，阶段切换中）'
            if event.get('event')=='finished' and event.get('returncode',0)!=0:current+='；返回值='+str(event['returncode'])
        f=Path(event.get('log',''))
        log=str(f) if f.is_file() else None
        if log:modified=datetime.fromtimestamp(f.stat().st_mtime).isoformat(timespec='seconds')
    return {'id':candidate['id'],'status':status,'gfn2':gfn2,'dft_inputs':manifest.get('dft_selected',len(inputs)) if manifest else len(inputs),
            'optimized':len(stages.get('opt',set())),'positive_frequencies':len(stages.get('freq',set())),
            'rejected_frequencies':len(rejected),'single_points':len(stages.get('sp',set())),
            'nmr':nmr,'nmr_tasks':sum(nmr.values()),'nmr_target_tasks':target,
            'nmr_complete_conformers':len(completed),'coverage_percent':100*coverage if coverage is not None else None,
            'current':current,'last_event_time':event.get('time') if event else None,
            'log':log,'log_modified':modified,'warnings':warnings}


def snapshot(project,cfg=None,output=None):
    project=Path(project).resolve()
    if cfg is None:
        f=next((project/n for n in ('job.json','job_compute.json') if (project/n).is_file()),None)
        if not f:raise ValueError('找不到项目计算配置')
        cfg=json.loads(f.read_text(encoding='utf-8'))
    output=Path(output) if output else (project/cfg['output']).resolve()
    events=recent_events(project)
    return {'project':str(project),'checked_at':datetime.now().isoformat(timespec='seconds'),
            'candidates':[candidate_snapshot(c,cfg,output,events) for c in cfg['candidates']]}


def format_snapshot(data):
    lines=['项目：'+data['project'],'读取时间：'+data['checked_at']]
    for r in data['candidates']:
        coverage='{:.2f}%'.format(r['coverage_percent']) if r['coverage_percent'] is not None else '待最终汇总'
        tasks=str(r['nmr_tasks'])+'/'+str(r['nmr_target_tasks']) if r['nmr_target_tasks'] is not None else str(r['nmr_tasks'])+'（目标待汇总）'
        lines.append('结构 {}：{}；GFN2={}，DFT输入={}，优化完成={}，正频率={}，ORCA完成={}，覆盖={}'.format(
            r['id'],r['status'],r['gfn2'] if r['gfn2'] is not None else '待生成',r['dft_inputs'],r['optimized'],r['positive_frequencies'],r['single_points'],coverage))
        lines.append('  NMR任务已完成：{}；{}；三方法全部完成的构象={}。'.format(
            tasks,'，'.join('NMR{}={}'.format(m,n) for m,n in r['nmr'].items()),r['nmr_complete_conformers']))
        lines.append('  当前：'+r['current'])
        if r['log']:lines.append('  独立日志：{}；最后更新：{}'.format(r['log'],r['log_modified']))
        if r['rejected_frequencies']:lines.append('  频率剔除构象数：'+str(r['rejected_frequencies']))
        lines.extend('  提示：'+w for w in r['warnings'])
    lines.append('数字按已完成回执计数；运行中的单点尚不计入。覆盖率仅针对已计算集合；进程是否仍运行请结合队列。')
    return '\n'.join(lines)


def queue_text(project):
    warnings=[];history=read_json(Path(project)/'提交记录.json',warnings)
    submission=read_json(Path(project)/'SUBMISSION.json',warnings)
    job=history[-1].get('job_id') if isinstance(history,list) and history else (submission or {}).get('job_id')
    if not job:return '尚无提交记录；未查询或提交计算。'
    if submission and submission.get('job_id')!=job:return '提交记录与SUBMISSION作业号不一致，请核查。'
    try:
        r=subprocess.run(['qstat','-f',str(job)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True,timeout=15)
    except (OSError,subprocess.TimeoutExpired) as exc:return '队列查询暂不可用：'+str(exc)
    if r.returncode:
        if re.search(r'Unknown Job|Unknown job|Invalid Job',r.stdout):return '作业 {} 已不在调度器记录中；请结合PBS结束日志核实退出状态。'.format(job)
        return '队列查询失败：'+r.stdout.strip()
    keep=[line.strip() for line in r.stdout.splitlines() if any(k in line for k in ('job_state =','resources_used.walltime =','resources_used.cput =','exit_status ='))]
    return '队列作业：'+str(job)+'\n'+'\n'.join(keep)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('project')
    parser.add_argument('--html',action='store_true',help='只刷新项目内的实时进度HTML快照')
    args=parser.parse_args(argv)
    data=snapshot(args.project)
    print(format_snapshot(data))
    print(queue_text(args.project))
    if args.html:
        project=Path(data['project'])
        text='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>实时计算进度</title><style>body{font:16px sans-serif;margin:32px;line-height:1.6}pre{white-space:pre-wrap;overflow-wrap:anywhere}</style><h1>实时计算进度</h1><p>本页是02入口最近一次读取的快照；再次运行02后刷新本页。完整文件索引见计算总览.html。</p><pre>'+html.escape(format_snapshot(data))+'</pre></html>'
        fd,name=tempfile.mkstemp(prefix='.live_progress_',dir=str(project))
        try:
            with os.fdopen(fd,'w',encoding='utf-8') as f:f.write(text)
            os.replace(name,project/'实时进度.html')
        finally:
            if os.path.exists(name):os.unlink(name)
        print('已刷新 实时进度.html；此页是本次读取快照。')
    return 0


if __name__=='__main__':
    try:raise SystemExit(main())
    except (OSError,ValueError,KeyError) as exc:raise SystemExit('进度读取失败（没有提交或修改计算）：'+str(exc))
