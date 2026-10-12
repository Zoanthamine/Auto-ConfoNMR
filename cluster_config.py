"""Site-specific PBS environment, separate from the scientific protocol."""
import argparse
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import sys

ROOT = Path(__file__).resolve().parent
EXTERNAL_RESOURCES = (
    'isostat', 'sobtop', 'sobtop.ini', 'atomtype', 'mol2_xyz', 'xyz2QC',
    'ATOMTYPE_AMBER.DEF', 'ATOMTYPE_GFF.DEF', 'CORR_NAME_TYPE.DAT',
    'LJ_param.dat', 'assign_AT.dat', 'bondcrit.dat', 'bonded_param.dat',
)
SOFTWARE = ('g16', 'formchk', 'orca', 'openmpi', 'crest', 'shermo',
            'multiwfn', 'obabel', 'vmd', 'gmx')
EXTERNAL_WORKBOOKS = ('DP4plus.xlsx', 'STS.xlsx')


def missing_resources(site):
    source = Path(site['resource_bundle'])
    names = list(EXTERNAL_RESOURCES) + ['upstream/'+n for n in EXTERNAL_WORKBOOKS]
    return [name for name in names if not (source/name).is_file()]


def load_site(path=None):
    """Read explicit local settings; never silently use the author's cluster."""
    path = Path(path or ROOT/'config/cluster.json')
    if not path.is_file():
        raise ValueError('缺少集群配置：复制 config/cluster.example.json 为 config/cluster.json，按迁移说明填写。')
    site = json.loads(path.read_text(encoding='utf-8-sig'))
    required = {'python', 'projects_root', 'scratch_root', 'resource_bundle',
                'software', 'modules', 'gaussian_profile', 'path_prepend',
                'ld_library_path_prepend', 'scheduler'}
    if set(site) != required:
        raise ValueError('集群配置字段有缺失或多余，请对照 cluster.example.json')
    if set(site['software']) != set(SOFTWARE):
        raise ValueError('software 字段必须包含：'+', '.join(SOFTWARE))
    for key in ['python', 'projects_root', 'scratch_root', 'resource_bundle']:
        site[key] = str(Path(site[key]).expanduser())
        if not Path(site[key]).is_absolute() or 'CHANGE_ME' in site[key]:
            raise ValueError('请填写绝对路径：'+key)
    for key, value in site['software'].items():
        if not isinstance(value, str) or not value or 'CHANGE_ME' in value:
            raise ValueError('请填写软件路径或可执行命令：'+key)
        site['software'][key] = os.path.expanduser(value)
    for key in ['modules', 'path_prepend', 'ld_library_path_prepend']:
        if not isinstance(site[key], list) or not all(isinstance(v, str) and '\n' not in v for v in site[key]):
            raise ValueError('应为不含换行的字符串列表：'+key)
    scheduler = site['scheduler']
    if set(scheduler) != {'type', 'queue', 'walltime', 'node_resource', 'extra_qsub_args'} or scheduler['type'] != 'pbs':
        raise ValueError('本版提供PBS提交器；Slurm/LSF需要独立适配，不能直接填写为pbs')
    if (not scheduler['queue'] or 'CHANGE_ME' in scheduler['queue']
            or not re.fullmatch(r'\d+:\d{2}:\d{2}', scheduler['walltime'])
            or '{cores}' not in scheduler['node_resource']):
        raise ValueError('请填写PBS队列、时限和包含{cores}的节点资源格式')
    if any('\n' in str(v) or '\r' in str(v) for v in [scheduler['queue'], scheduler['walltime'], scheduler['node_resource']]):
        raise ValueError('PBS字段不能含换行')
    if site['gaussian_profile'] and (not Path(site['gaussian_profile']).expanduser().is_absolute()
                                    or 'CHANGE_ME' in site['gaussian_profile']):
        raise ValueError('gaussian_profile应为空或真实绝对路径')
    if not Path(site['software']['openmpi']).is_absolute():
        raise ValueError('openmpi必须是安装根目录的绝对路径')
    if not isinstance(scheduler['extra_qsub_args'], list) or not all(isinstance(v, str) and '\n' not in v for v in scheduler['extra_qsub_args']):
        raise ValueError('extra_qsub_args应为参数列表')
    return site


def install_resources(site, workflow):
    """Copy licensed, locally installed helper resources into the fixed snapshot."""
    source = Path(site['resource_bundle'])
    missing = missing_resources(site)
    if missing:
        raise ValueError('外部资源包缺少：'+', '.join(missing)+'；参见 README.md 的迁移配置部分')
    target = Path(workflow)/'resources'
    target.mkdir(parents=True, exist_ok=True)
    for name in EXTERNAL_RESOURCES:
        shutil.copy2(source/name, target/name)
    upstream = Path(workflow)/'upstream'
    upstream.mkdir(exist_ok=True)
    for name in EXTERNAL_WORKBOOKS:
        shutil.copy2(source/'upstream'/name, upstream/name)


def runtime_text(site):
    """Generate quoted shell statements from trusted local configuration."""
    q = shlex.quote
    lines = ['#!/usr/bin/env bash', '# Generated from config/cluster.json; no scientific settings are changed.',
             'export CONFONMR_PYTHON='+q(site['python'])]
    if site['modules']:
        lines += ['if ! type module >/dev/null 2>&1; then',
                  "  printf 'module命令不可用，请按迁移说明初始化Environment Modules/Lmod。\\n' >&2", '  return 1', 'fi']
        lines += ['module load '+q(name) for name in site['modules']]
    if site['path_prepend']:
        lines += ['export PATH='+q(':'.join(site['path_prepend']))+':"$PATH"']
    if site['ld_library_path_prepend']:
        lines += ['export LD_LIBRARY_PATH='+q(':'.join(site['ld_library_path_prepend']))+':"${LD_LIBRARY_PATH:-}"']
    if site['gaussian_profile']:
        profile = str(Path(site['gaussian_profile']).expanduser())
        # Gaussian's vendor profile can inspect unset variables or return an
        # architecture probe code; conditionally source it, then validate output.
        lines += ['export g16root='+q(str(Path(site['software']['g16']).parent.parent)),
                  'set +u', '# shellcheck disable=SC1090,SC1091', 'if source '+q(profile)+'; then :; else',
                  "  printf 'Gaussian环境初始化失败。\\n' >&2", '  return 1', 'fi', 'set -u',
                  ': "${GAUSS_EXEDIR:?Gaussian profile did not set GAUSS_EXEDIR}"']
    lines += ['export Multiwfnpath='+q(str(Path(site['software']['multiwfn']).parent)), '']
    return '\n'.join(lines)


def scheduler_arguments(site, cfg):
    s = site['scheduler']
    resource = s['node_resource'].format(cores=cfg['cores'], memory_gb=cfg['memory_gb'])
    args = ['-q', s['queue'], '-l', resource, '-l', 'walltime='+s['walltime'], *s['extra_qsub_args']]
    if any('\n' in a or '\r' in a for a in args):
        raise ValueError('PBS参数不能含换行')
    return args


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['check', 'scheduler-args'])
    p.add_argument('job', nargs='?')
    a = p.parse_args()
    site = load_site()
    if a.action == 'scheduler-args':
        if not a.job:p.error('需要 job.json')
        cfg = json.loads(Path(a.job).read_text())
        print('\n'.join(scheduler_arguments(site, cfg)))
    else:
        missing = missing_resources(site)
        if missing:raise ValueError('外部资源缺少：'+', '.join(missing))
        if not Path(site['python']).is_file():raise ValueError('python路径不存在')
        print('集群配置格式与资源包完整性通过；软件运行环境仍需在目标计算节点验证。')


if __name__ == '__main__':
    try:main()
    except (ValueError, OSError, KeyError) as exc:sys.exit('停止：'+str(exc))
