"""Versioned choices for new projects; old projects keep their original identity."""
import copy
import math
import re

LAYOUT='dp4_dp4_sts_v4'
NMR_PROFILES={
    'legacy_pcm_g09':('IEFPCM','G09Defaults','IEFPCM + G09Defaults'),
    'pcm_g16':('IEFPCM','G16','PCM/IEFPCM + G16默认值（版本差异对照）'),
    'smd_g16':('SMD','G16','SMD + G16默认值（探索对照）'),
}
# Edis: Hartree energies are read by isostat; its -Edis option uses kcal/mol.
# Gdis: RMSD/Angstrom. T: isostat population display temperature, not Shermo T.
CLUSTERS={
    'coarse':{
        'GFN1':{'Edis':.5,'Gdis':.5,'T':300,'Nout':100,'Eout':100},
        'GFN2':{'Edis':.3,'Gdis':.2,'T':300,'Nout':100,'Eout':100},
        'DFT':{'Edis':.3,'Gdis':.2,'T':300,'Nout':100,'Eout':100},
    },
    'fine':{
        'GFN1':{'Edis':.2,'Gdis':.2,'T':350,'Nout':100,'Eout':100},
        'GFN2':{'Edis':.2,'Gdis':.1,'T':350,'Nout':100,'Eout':100},
        'DFT':{'Edis':.2,'Gdis':.1,'T':350,'Nout':100,'Eout':100},
    },
}


def method_map(cfg):
    if cfg.get('method_layout')==LAYOUT:
        return {'01':('B3LYP','6-31G(d)','dp4plus'),
                '02':('mPW1PW91','6-311+G(d,p)','dp4plus'),
                '03':('wB97XD','6-31G(d)','sts')}
    if cfg.get('method_layout') is not None:raise ValueError('未知方法编号版本')
    return {'01':('B3LYP','6-31G(d)','dp4plus'),
            '02':('wB97XD','6-31G(d)','sts'),
            '03':('mPW1PW91','6-311+G(d,p)','dp4plus')}


def template_method(method,cfg):
    # Preserve the method order recorded by configurations without a layout key.
    method_map(cfg)
    return method if cfg.get('method_layout')==LAYOUT else {'01':'01','02':'03','03':'02'}[method]


def apply_choices(cfg,nmr_profile='legacy_pcm_g09',clustering_profile='coarse'):
    if nmr_profile not in NMR_PROFILES:raise ValueError('未知NMR方案：'+nmr_profile)
    if clustering_profile not in CLUSTERS:raise ValueError('未知聚类档位：'+clustering_profile)
    cfg['protocol']='Auto_ConfoNMR_V1.1_custom'
    cfg['method_layout']=LAYOUT
    cfg['nmr_profile']=nmr_profile
    model,defaults,_=NMR_PROFILES[nmr_profile]
    cfg['nmr_settings']={'solvent_model':model,'gaussian_defaults':defaults}
    cfg['clustering_profile']=clustering_profile
    cfg['clustering']=copy.deepcopy(CLUSTERS[clustering_profile])
    common=['Geometry uses gas B3LYP-D3BJ/6-31G(d); not original DP4+ B3LYP/6-31G(d) or STS solvent B3LYP-D3BJ/TZVP',
            'Weights use gas ORCA energies + Shermo quasi-RRHO; DP4+ originally uses corresponding NMR single-point energies at 298 K; STS uses solvent B3LYP-D3BJ/TZVP Gibbs energies',
            'Fixed DFT budget and legacy Sobtop guessed torsions limit conformer coverage']
    if model=='SMD':common.append('NMR uses SMD chloroform; calibration tables require PCM/IEFPCM')
    if defaults=='G16':common.append('Gaussian 16 defaults differ from Gaussian 09 reference calculations')
    cfg['postprocessing']=[
        {'kind':kind,'method':mid,'parameter_key':functional+'|PCM|'+basis if kind=='dp4plus' else 'CDCl3',
         'allow_exploratory':True,'differences':list(common)}
        for mid,(functional,basis,kind) in method_map(cfg).items()]
    return cfg


def nmr_text(text,cfg):
    settings=cfg.get('nmr_settings')
    if settings is None:
        # A configuration without explicit settings retains its SMD/G16 semantics.
        if re.search(r'(?i)\bscrf\s*=?\s*\(smd\b',text) and not re.search(r'(?i)\bg09defaults\b',text):return text
        settings={'solvent_model':'SMD','gaussian_defaults':'G16'}
    model=settings['solvent_model'];defaults=settings['gaussian_defaults']
    if model not in ('IEFPCM','SMD') or defaults not in ('G09Defaults','G16'):
        raise ValueError('不支持的NMR溶剂模型/兼容设置')
    pattern=r'(?i)\bscrf\s*=?\s*\([^\n)]*\)'
    if len(re.findall(pattern,text))!=1:raise ValueError('NMR模板应只有一个SCRF配置')
    text=re.sub(pattern,'SCRF=('+model+',Solvent=Chloroform)',text)
    text=re.sub(r'(?i)\s+g09defaults\b','',text)
    if defaults=='G09Defaults':
        text=re.sub(r'(?m)^(#[^\n]*)$',r'\1 G09Defaults',text)
    return text


def cluster_options(cfg,stage):
    c=cfg['clustering']
    if stage in c:
        values=c[stage]
        return [part for key in ('Edis','Gdis','T','Nout','Eout') for part in ('-'+key,values[key])]
    # Preserve historic callers exactly; their frozen projects used T350.
    return ['-Edis',c['Edis'],'-Gdis',c['GFN1_Gdis'] if stage=='GFN1' else c['Gdis'],'-T',350]


def validate_choices(cfg):
    method_map(cfg)
    if cfg.get('method_layout')==LAYOUT and not {'nmr_settings','nmr_profile','clustering_profile'} <= set(cfg):
        raise ValueError('配置缺少显式NMR或聚类方案，拒绝回退到历史默认值')
    if 'nmr_settings' in cfg:
        profile=cfg.get('nmr_profile')
        if profile not in NMR_PROFILES:raise ValueError('未知NMR方案')
        model,defaults,_=NMR_PROFILES[profile]
        if cfg['nmr_settings']!={'solvent_model':model,'gaussian_defaults':defaults}:
            raise ValueError('NMR方案与记录的实际参数不一致')
    c=cfg['clustering']
    if set(c)=={'GFN1','GFN2','DFT'}:
        for stage in c:
            v=c[stage]
            if set(v)!={'Edis','Gdis','T','Nout','Eout'}:raise ValueError('聚类字段遗漏')
            if not all(isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x) and x>0 for x in v.values()):
                raise ValueError('聚类阈值必须为有限正数')
            if not isinstance(v['Nout'],int):raise ValueError('Nout必须为整数')
        name=cfg.get('clustering_profile')
        if name not in CLUSTERS or c!=CLUSTERS[name]:raise ValueError('聚类档位与实际阈值不一致')
    elif not all(isinstance(v,(int,float)) and math.isfinite(v) and v>0 for v in c.values()):
        raise ValueError('旧聚类阈值无效')


def overview(cfg):
    lines=['# 本项目计算方案','',
           'NMR方案：'+NMR_PROFILES.get(cfg.get('nmr_profile'),('','','历史SMD方案'))[2],
           '聚类档位：'+cfg.get('clustering_profile','历史固定阈值'),'',
           '|编号|NMR泛函/基组|后处理|','|---|---|---|']
    for mid,(f,b,k) in method_map(cfg).items():lines.append('|{}|{}/{}|{}|'.format(mid,f,b,'DP4+' if k=='dp4plus' else 'STS（仅13C）'))
    lines+=['','## 聚类参数','', '|阶段|Edis/kcal mol-1|Gdis/Å|isostat T/K|Nout|Eout/kcal mol-1|','|---|---:|---:|---:|---:|---:|']
    for stage,v in cfg['clustering'].items():
        if isinstance(v,dict):lines.append('|{}|{}|{}|{}|{}|{}|'.format(stage,v['Edis'],v['Gdis'],v['T'],v['Nout'],v['Eout']))
    lines+=['','isostat的T只影响该工具分布显示；最终混合比例使用Shermo温度{} K。'.format(cfg['shermo']['T']),
            '最大DFT构象预算{}；累计分布{}仅针对完成的DFT集合。'.format(cfg['max_conformers'],cfg['population_coverage']),
            '共同几何为气相B3LYP-D3BJ/6-31G(d)，权重来自气相ORCA/Shermo。共同几何与权重条件不同于方法论文，须保留适用性差异记录。',
            'DP4+/STS输出保留自定义方案标记。详细差异：']
    lines+=['- '+s for s in cfg['postprocessing'][0]['differences']]
    return '\n'.join(lines)+'\n'
