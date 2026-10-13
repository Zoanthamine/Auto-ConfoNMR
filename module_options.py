"""Explicit V2 module choices. Old configurations retain the V1 calculation path."""
import copy
import math

DEFAULT_MODULES = {'base':'nmr', 'conformer_weight':'free_energy', 'coupling':False, 'ecd':False,
                   'ecd_methods':['01','02','03'], 'ecd_states':30, 'ecd_solvent':'methanol'}
DEFAULT_ANALYSIS = {'top_n':6, 'coupling_allow_exploratory':False,
                    'j_experiment':'experiment_j.csv', 'j_assignments':'assignments_j.csv',
                    'ecd_experiment':'experiment_ecd.csv', 'ecd_fwhm_ev':0.5,
                    'ecd_gauge':'velocity', 'ecd_energy_scale':1.0, 'ecd_shift_nm':0.0,
                    'ecd_min_nm':150.0, 'ecd_max_nm':800.0, 'ecd_step_nm':1.0,
                    'ecd_normalize':False}


def options(cfg):
    result=copy.deepcopy(DEFAULT_MODULES)
    result.update(cfg.get('modules',{}))
    return result


def analysis(cfg):
    result=copy.deepcopy(DEFAULT_ANALYSIS)
    result.update(cfg.get('module_analysis',{}))
    return result


def validate(cfg):
    from nmrflow import require
    m,a=options(cfg),analysis(cfg)
    require(set(m)==set(DEFAULT_MODULES) and set(a)==set(DEFAULT_ANALYSIS),'Unknown V2 module setting')
    require(m['base'] in ('nmr','conformers'),'Base must be nmr or conformers')
    require(m['conformer_weight'] in ('gfn2_energy','free_energy'),'Unknown conformer population model')
    require(all(isinstance(m[k],bool) for k in ('coupling','ecd')),'Module switches must be boolean')
    require(m['conformer_weight']=='free_energy' or (m['base']=='conformers' and not m['coupling'] and not m['ecd']),
            'GFN2 electronic-energy mode is only for standalone conformers')
    require(m['ecd_solvent'] in ('methanol','chloroform','acetonitrile','water','DMSO'), 'Unsupported explicit ECD solvent')
    require(m['ecd_methods'] and set(m['ecd_methods']) <= {'01','02','03'} and len(set(m['ecd_methods']))==len(m['ecd_methods']), 'Invalid ECD methods')
    require(isinstance(m['ecd_states'],int) and not isinstance(m['ecd_states'],bool) and m['ecd_states']>0,'Invalid ECD state count')
    require(not m['ecd'] or cfg['multiplicity']==1,'This ECD module currently supports closed-shell singlet TDDFT only')
    require(isinstance(a['top_n'],int) and not isinstance(a['top_n'],bool) and a['top_n']>0,'Export count must be a positive integer')
    require(isinstance(a['coupling_allow_exploratory'],bool) and isinstance(a['ecd_normalize'],bool),'Invalid analysis switch')
    require(a['ecd_gauge'] in ('length','velocity'),'Invalid ECD rotatory-strength gauge')
    require(all(math.isfinite(a[k]) for k in ('ecd_fwhm_ev','ecd_energy_scale','ecd_shift_nm','ecd_min_nm','ecd_max_nm','ecd_step_nm')),'Nonfinite ECD plotting value')
    require(a['ecd_fwhm_ev']>0 and a['ecd_energy_scale']>0 and 0<a['ecd_min_nm']<a['ecd_max_nm'] and a['ecd_step_nm']>0,'Invalid ECD plotting range/broadening')
    require((a['ecd_max_nm']-a['ecd_min_nm'])/a['ecd_step_nm']<=100000,'ECD grid is unnecessarily large')


def ecd_text(text,cfg):
    import re
    m=options(cfg)
    return re.sub(r'nstates=\d+', 'nstates='+str(m['ecd_states']),
                  re.sub(r'solvent=methanol', 'solvent='+m['ecd_solvent'],text,flags=re.I),flags=re.I)


def settings_text(cfg):
    m,a=options(cfg),analysis(cfg)
    return '''
[模块]
# nmr=现有三分支核磁；conformers=只算构象。J/ECD可附加到精修模式。
base = {base}
# free_energy=DFT优化/频率+ORCA+Shermo；gfn2_energy=快速GFN2电子能估算。
# 快速档只用于单独构象，不能与J或ECD混用。
conformer_weight = {conformer_weight}
coupling = {coupling}
ecd = {ecd}
# ECD独立于核磁溶剂设置。默认保留旧三套模板，SMD甲醇，30个激发态。
ecd_methods = {ecd_methods}
ecd_states = {ecd_states}
ecd_solvent = {ecd_solvent}

[模块后处理]
# 默认导出权重最高的6个构象。只影响导出，不会追加DFT计算。
top_n = {top_n}
# 自定义DFT几何/Shermo布居不匹配dJ原文MMFF校准；默认只输出表。
# 审核后改true，可在有完整C/H/J实验归属时做明确标注的探索比较。
coupling_allow_exploratory = {coupling_allow_exploratory}
j_experiment = {j_experiment}
j_assignments = {j_assignments}
ecd_experiment = {ecd_experiment}
# 高斯展宽半高全宽FWHM，单位eV，不是sigma标准差。
ecd_fwhm_ev = {ecd_fwhm_ev}
ecd_gauge = {ecd_gauge}
# 下列是显式显示修正；默认不校正。修改会保留原始曲线。
ecd_energy_scale = {ecd_energy_scale}
ecd_shift_nm = {ecd_shift_nm}
ecd_min_nm = {ecd_min_nm}
ecd_max_nm = {ecd_max_nm}
ecd_step_nm = {ecd_step_nm}
# false保留Δε物理单位；true额外生成归一化图，纵轴为任意单位。
ecd_normalize = {ecd_normalize}
'''.format(**dict(m,ecd_methods=','.join(m['ecd_methods'])),**a)


def read_settings(parser,cfg):
    if '模块' not in parser:
        return
    if set(parser['模块'])!=set(DEFAULT_MODULES) or set(parser['模块后处理'])!=set(DEFAULT_ANALYSIS):
        raise ValueError('模块配置字段有遗漏/拼写错误')
    m={k:parser.get('模块',k).strip() for k in DEFAULT_MODULES}
    for k in ('coupling','ecd'):m[k]=parser.getboolean('模块',k)
    m['ecd_methods']=[v.strip() for v in m['ecd_methods'].split(',')]
    m['ecd_states']=parser.getint('模块','ecd_states')
    a={k:parser.get('模块后处理',k).strip() for k in DEFAULT_ANALYSIS}
    for k in ('coupling_allow_exploratory','ecd_normalize'):a[k]=parser.getboolean('模块后处理',k)
    a['top_n']=parser.getint('模块后处理','top_n')
    for k in ('ecd_fwhm_ev','ecd_energy_scale','ecd_shift_nm','ecd_min_nm','ecd_max_nm','ecd_step_nm'):a[k]=parser.getfloat('模块后处理',k)
    cfg['modules'],cfg['module_analysis']=m,a
    validate(cfg)
