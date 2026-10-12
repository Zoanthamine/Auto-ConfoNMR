#!/usr/bin/env python3
"""Stage-isolated Huangjt NMR workflow. Python >=3.8; numpy required."""
import argparse
import uuid
import socket
import errno
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
import tempfile
from contextlib import contextmanager
from protocol_options import cluster_options, method_map, nmr_text, template_method, validate_choices

ROOT = Path(__file__).resolve().parent
ELEMENTS = ('X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn '
            'Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd '
            'In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu '
            'Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn').split()
NUMBER = r'[-+]?\d+(?:\.\d*)?(?:[DEde][-+]?\d+)?'


def require(ok, message):
    if not ok:
        raise ValueError(message)


@contextmanager
def output_lock(output):
    """One atomic shared gate on every node; kernel locks provide an extra check."""
    output.mkdir(parents=True, exist_ok=True)
    guard=output/'.workflow.guard'
    try:
        guard.mkdir()
    except FileExistsError as exc:
        raise ValueError('Another workflow/postprocessor is using '+str(output)+'; inspect '+str(guard/'owner.json')) from exc
    owner_path=guard/'owner.json'
    owner={'token':uuid.uuid4().hex,'pid':os.getpid(),'host':socket.gethostname(),
           'pbs_job_id':os.environ.get('PBS_JOBID'),'started':time.strftime('%Y-%m-%d %H:%M:%S'),
           'mechanism':'atomic_directory_and_kernel','kernel_lock_errno':None}
    try:
        save(owner_path,owner)
        with (output/'.workflow.lock').open('a+b') as handle:
            handle.seek(0)
            if os.name == 'nt':
                import msvcrt
                if os.fstat(handle.fileno()).st_size == 0:
                    handle.write(b'0');handle.flush()
                handle.seek(0)
                lock,unlock=lambda:msvcrt.locking(handle.fileno(),msvcrt.LK_NBLCK,1),lambda:msvcrt.locking(handle.fileno(),msvcrt.LK_UNLCK,1)
            else:
                import fcntl
                lock,unlock=lambda:fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB),lambda:fcntl.flock(handle,fcntl.LOCK_UN)
            acquired=False
            try:
                try:
                    lock();acquired=True
                except OSError as exc:
                    if os.name!='nt' and exc.errno in (errno.ENOLCK,errno.ENOSYS,errno.EOPNOTSUPP):
                        # All clients of this lock implementation still hold the SAME atomic-directory gate.
                        # Never fall back to a separate per-node/local lock.
                        owner.update(mechanism='atomic_directory',kernel_lock_errno=exc.errno)
                        save(owner_path,owner)
                    elif exc.errno in (errno.EACCES,errno.EAGAIN):
                        raise ValueError('Another workflow/postprocessor is using '+str(output)) from exc
                    else:
                        raise ValueError('Workflow lock service failed for '+str(output)+': '+str(exc)) from exc
                yield
            finally:
                if acquired:unlock()
    finally:
        if owner_path.is_file():
            require(json.loads(owner_path.read_text())['token']==owner['token'],'Workflow gate ownership changed; retained '+str(guard))
            owner_path.unlink()
        guard.rmdir()


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def calculation_config(cfg):
    return {k:v for k,v in cfg.items() if k not in ('postprocessing','experiment','assignments')}


def save(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    tmp.replace(path)


def xyz_read(path):
    lines = Path(path).read_text().splitlines()
    result = []
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i])
        require(n > 1 and i+n+2 <= len(lines), 'Invalid/truncated XYZ')
        geom = []
        for line in lines[i+2:i+n+2]:
            p = line.split()
            require(len(p) == 4 and p[0] in ELEMENTS[1:], 'Invalid XYZ atom')
            xyz = [float(x) for x in p[1:]]
            require(all(math.isfinite(x) for x in xyz), 'Nonfinite coordinate')
            geom.append([p[0]] + xyz)
        result.append((lines[i+1], geom))
        i += n+2
    require(bool(result), 'Empty XYZ')
    symbols = [x[0] for x in result[0][1]]
    require(all([x[0] for x in g] == symbols for _, g in result), 'Atom order differs across XYZ blocks')
    return result


def xyz_text(geom, comment=''):
    return '{}\n{}\n{}\n'.format(len(geom), comment, '\n'.join(
        '{} {:.10f} {:.10f} {:.10f}'.format(*x) for x in geom))


def ensemble_energy(comment):
    """Read CREST's Hartree-first or isostat's explicitly labelled comment."""
    labelled = re.fullmatch(r'\s*Energy\s*=\s*('+NUMBER+r')\s+(?:a\.u\.|Eh|Hartree)'
                            r'(?:\s+#Cluster:\s*\d+)?\s*',comment,re.I)
    if labelled:
        value = labelled[1]
    else:
        require(not re.search(r'\b(?:kcal|kJ|eV)\b',comment,re.I),
                'Ensemble energy must be in Hartree, not another unit')
        tokens = comment.split()
        require(tokens and re.fullmatch(NUMBER,tokens[0]),
                'Unsupported ensemble energy comment; expected Hartree-first CREST or Energy = ... a.u. isostat format')
        value = tokens[0]
    energy = float(value.replace('D','E').replace('d','e'))
    require(math.isfinite(energy),'Nonfinite ensemble energy')
    return energy


def select_ensemble(blocks, limit):
    indexed = [(i,c,g) for i,(c,g) in enumerate(blocks,1)]
    if limit and len(indexed) > limit:
        energies = {i:ensemble_energy(c) for i,c,_ in indexed}
        indexed.sort(key=lambda r:(energies[r[0]],r[0]))
        indexed = indexed[:limit]
    return indexed


def geometry_text(geom):
    return '\n'.join('{} {:.10f} {:.10f} {:.10f}'.format(*a) for a in geom)


def gaussian_read(path, stage):
    text = Path(path).read_text(errors='replace')
    require('Normal termination of Gaussian' in text and 'Error termination' not in text,
            '{}: Gaussian did not terminate normally'.format(path))
    energies = re.findall(r'SCF Done:.*?=\s*(' + NUMBER + ')', text)
    require(bool(energies), '{}: missing SCF energy'.format(path))
    orientations = re.findall(r'(?:Standard|Input) orientation:\s*\n\s*-+\n.*?\n\s*-+\n(.*?)\n\s*-+', text, re.S)
    require(bool(orientations), '{}: missing atomic geometry'.format(path))
    geom = []
    for row in orientations[-1].splitlines():
        p = row.split()
        require(len(p) == 6 and 0 < int(p[1]) < len(ELEMENTS), 'Unsupported orientation table')
        geom.append([ELEMENTS[int(p[1])]] + [float(v) for v in p[3:6]])
    frequencies = [float(x.replace('D', 'E')) for row in re.findall(r'Frequencies --\s*([^\n]+)', text)
                   for x in row.split()]
    if stage == 'opt':
        require('Optimization completed.' in text, '{}: optimization not converged'.format(path))
    if stage == 'freq':
        require(bool(frequencies), '{}: no frequencies'.format(path))
        require(len(frequencies) in (3*len(geom)-6,3*len(geom)-5), '{}: incomplete frequency list'.format(path))
        require(all(math.isfinite(v) for v in frequencies), 'Nonfinite frequency')
        require(min(frequencies) > 0, '{}: imaginary/zero frequencies {}'.format(path, min(frequencies)))
    shielding = {}
    for a, el, val in re.findall(r'^\s*(\d+)\s+([A-Z][a-z]?)\s+Isotropic\s*=\s*(' + NUMBER + ')', text, re.M):
        idx = int(a)
        require(idx not in shielding, '{}: multiple shielding tables; split Link1 jobs first'.format(path))
        shielding[idx] = (el, float(val.replace('D', 'E')))
    if stage == 'nmr':
        require(len(shielding) == len(geom), '{}: incomplete shielding table'.format(path))
        require(all(shielding[i+1][0] == a[0] for i, a in enumerate(geom)), 'Shielding atom identity mismatch')
    require(all(math.isfinite(x) for a in geom for x in a[1:]), 'Invalid Gaussian geometry')
    require(all(math.isfinite(v) for _,v in shielding.values()), 'Nonfinite shielding')
    require(math.isfinite(float(energies[-1].replace('D','E'))), 'Nonfinite SCF energy')
    return {'geometry': geom, 'energy': float(energies[-1].replace('D', 'E')),
            'frequencies': frequencies, 'shielding': shielding}


def rms(a, b):
    import numpy as np
    require([x[0] for x in a] == [x[0] for x in b], 'Atom identity/order mismatch')
    p, q = np.array([x[1:] for x in a],dtype=float), np.array([x[1:] for x in b],dtype=float)
    p -= p.mean(axis=0)
    q -= q.mean(axis=0)
    u, _, vt = np.linalg.svd(p.T @ q)
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(u @ vt)
    return float(np.sqrt(((p @ (u @ correction @ vt) - q)**2).sum()/len(a)))


def command(argv, cwd, log, env=None, stdin=None):
    cwd = Path(cwd)
    cwd.mkdir(parents=True, exist_ok=True)
    stream = os.environ.get('NMR_STREAM_LOG') == '1'
    detail = stream and (env if env is not None else os.environ).get('NMR_DETAIL_LOG', '1') != '0'
    project = os.environ.get('NMR_PROJECT_DIR')
    if stream:
        print('\n=== 软件开始 {} | {} ===\n工作目录：{}\n独立日志：{}'.format(
            time.strftime('%Y-%m-%d %H:%M:%S'), ' '.join(str(a) for a in argv), cwd, log), flush=True)
    with Path(log).open('wb') as out:
        if project:
            from project_views import event, refresh_from_env
            event(project,'started',argv=[str(a) for a in argv],cwd=str(cwd),log=str(log))
            refresh_from_env()
        if detail:
            # A temporary input file avoids pipe deadlocks for interactive batch menus.
            with tempfile.TemporaryFile() as inp:
                if stdin is not None:
                    inp.write(stdin.encode()); inp.seek(0)
                with subprocess.Popen([str(a) for a in argv],cwd=str(cwd),env=env,stdin=inp,
                                      stdout=subprocess.PIPE,stderr=subprocess.STDOUT) as process:
                    while True:
                        chunk = process.stdout.read1(65536)
                        if not chunk:
                            break
                        out.write(chunk); out.flush()
                        if hasattr(sys.stdout,'buffer'):
                            sys.stdout.buffer.write(chunk); sys.stdout.buffer.flush()
                        else:
                            sys.stdout.write(chunk.decode('utf-8',errors='replace')); sys.stdout.flush()
                    code = process.wait()
        else:
            result = subprocess.run([str(a) for a in argv], cwd=str(cwd), env=env,
                                    input=stdin.encode() if stdin is not None else None,
                                    stdout=out, stderr=subprocess.STDOUT)
            code = result.returncode
    if project:
        event(project,'finished',argv=[str(a) for a in argv],cwd=str(cwd),log=str(log),returncode=code)
        refresh_from_env()
    if stream:
        print('\n=== 软件结束 {} | 返回值 {} | {} ==='.format(time.strftime('%Y-%m-%d %H:%M:%S'),code,log),flush=True)
    require(code == 0, 'Command failed ({}); see {}'.format(code, log))


def artifact_stem(directory, cfg, software, stage):
    """Human-readable filenames for new projects; frozen legacy names stay valid."""
    directory = Path(directory)
    if cfg.get('output_layout') != 'direct_v1' or not re.fullmatch(r'c\d+', directory.parent.name):
        return 'job'
    candidate, conformer = directory.parents[2].name, directory.parent.name.upper()
    require(re.fullmatch(r'[A-Za-z0-9_-]+',candidate), 'Invalid candidate artifact name')
    if stage == 'opt':
        attempt = re.fullmatch(r'opt([123])', directory.name)
        require(attempt, 'Missing optimization attempt number')
        return '{}_opt_{}_{}_A{:02d}'.format(software,candidate,conformer,int(attempt[1]))
    phase = directory.name if stage == 'nmr' else stage
    return '{}_{}_{}_{}'.format(software,phase,candidate,conformer)


def quantum(template, geom, cfg, directory, stage):
    directory.mkdir(parents=True, exist_ok=True)
    text = Path(template).read_text()
    if stage == 'nmr':
        text = nmr_text(text,cfg)
    text = re.sub(r'%mem=.*', '%mem={}GB'.format(cfg['memory_gb']), text, flags=re.I)
    text = re.sub(r'%nprocshared=.*', '%nprocshared={}'.format(cfg['cores']), text, flags=re.I)
    text = re.sub(r'^[ \t]*0[ \t]+1[ \t]*$', '{} {}'.format(cfg['charge'], cfg['multiplicity']), text, flags=re.M)
    text = text.replace('[GEOMETRY]', geometry_text(geom))
    stem = artifact_stem(directory, cfg, 'Gau', stage)
    text = '%chk='+stem+'.chk\n' + re.sub(r'^%chk=.*\n', '', text, flags=re.M | re.I)
    inp, log, receipt = directory/(stem+'.gjf'), directory/(stem+'.log'), directory/'receipt.json'
    signature = hashlib.sha256(text.encode()).hexdigest()
    if receipt.exists():
        old = json.loads(receipt.read_text())
        require(old['input_sha256'] == signature and digest(log) == old['output_sha256'],
                'Changed input/output in {}; choose a new run directory'.format(directory))
        return gaussian_read(log, stage)
    require(not log.exists(), 'Unfinished stage {}; archive it before retrying'.format(directory))
    inp.write_text(text)
    env = os.environ.copy()
    env['NMR_DETAIL_LOG'] = '0'
    scratch = directory/'scratch'
    scratch.mkdir(exist_ok=True)
    env['GAUSS_SCRDIR'] = str(scratch.resolve())
    command([cfg['software']['g16']], directory, log, env=env, stdin=text)
    try:
        result = gaussian_read(log, stage)
    except ValueError as exc:
        if stage == 'freq' and 'imaginary/zero frequencies' in str(exc):
            save(receipt, {'input_sha256':signature,'output_sha256':digest(log),'stage':stage,'status':'rejected_frequency'})
        raise
    require(rms(geom, result['geometry']) < (100 if stage == 'opt' else 1e-4), 'Single-point geometry changed')
    save(receipt, {'input_sha256': signature, 'output_sha256': digest(log), 'stage': stage})
    return result


def orca_sp(geom, cfg, directory):
    directory.mkdir(parents=True, exist_ok=True)
    text = (ROOT/'resources/05_tem.inp').read_text()
    text = re.sub(r'%maxcore\s+\d+', '%maxcore {}'.format(cfg['orca_maxcore_mb']), text, flags=re.I)
    text = re.sub(r'nprocs\s+\d+', 'nprocs {}'.format(cfg['cores']), text, flags=re.I)
    text = re.sub(r'\* xyz\s+0\s+1', '* xyz {} {}'.format(cfg['charge'], cfg['multiplicity']), text)
    text = text.replace('[GEOMETRY]', geometry_text(geom))
    stem = artifact_stem(directory, cfg, 'Orca', 'sp')
    inp, log = directory/(stem+'.inp'), directory/(stem+'.log')
    receipt = directory/'receipt.json'
    signature = hashlib.sha256(text.encode()).hexdigest()
    if receipt.exists():
        old = json.loads(receipt.read_text())
        require(old['input_sha256'] == signature and old['output_sha256'] == digest(log), 'ORCA stage changed')
    else:
        require(not log.exists(), 'Unfinished ORCA stage {}; archive before retry'.format(directory))
        inp.write_text(text)
        env = os.environ.copy()
        env['NMR_DETAIL_LOG'] = '0'
        mpi = cfg['software']['openmpi']
        orca_dir = str(Path(cfg['software']['orca']).parent)
        env['PATH'] = mpi+'/bin:'+orca_dir+':/usr/bin:/bin'
        env['LD_LIBRARY_PATH'] = mpi+'/lib:'+orca_dir+':'+env.get('LD_LIBRARY_PATH', '')
        env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', OMPI_MCA_btl='self,vader,tcp')
        command([cfg['software']['orca'], inp.name], directory, log, env=env)
    output = log.read_text(errors='replace')
    values = re.findall(r'FINAL SINGLE POINT ENERGY\s+('+NUMBER+')', output)
    require('ORCA TERMINATED NORMALLY' in output and len(values) == 1, 'ORCA output invalid')
    energy = float(values[0])
    require(math.isfinite(energy), 'Invalid ORCA energy')
    save(receipt, {'input_sha256': signature, 'output_sha256': digest(log)})
    return energy


def populations(records, cfg, directory):
    directory.mkdir(parents=True, exist_ok=True)
    sher = cfg['shermo']
    settings = (ROOT/'resources/shermo_settings.ini').read_text()
    allowed = {'T','P','ilowfreq','intpvib','ravib','sclZPE','sclheat','sclS','sclCV','conc','imode','defmass'}
    for key, val in sher.items():
        require(key in allowed, 'Unsupported Shermo setting: '+key)
        require('\n' not in str(val) and '\r' not in str(val), 'Invalid Shermo setting value')
        pattern = r'(?m)^[ \t]*'+re.escape(key)+r'[ \t]*=.*$'
        line = key+'= '+str(val)
        if re.search(pattern, settings):
            settings = re.sub(pattern, lambda _: line, settings)
        else:
            settings = line+'\n'+settings
    (directory/'settings.ini').write_text(settings)
    # Run individually: explicit ID-to-energy mapping, never parse sorted percentages.
    for row in records:
        name = 'Shermo_weight_{}_{}'.format(directory.parent.name,row['id'].upper()) if cfg.get('output_layout')=='direct_v1' else row['id']
        log = directory/(name+'.log')
        argv = [cfg['software']['shermo'], row['frequency'], '-E', str(row['sp_energy'])]
        command(argv, directory, log)
        matches = re.findall(r'Sum of electronic energy and thermal correction to G:\s*('+NUMBER+')', log.read_text())
        require(len(matches) == 1, 'Expected one Shermo Gibbs energy: '+str(log))
        row['gibbs_hartree'] = float(matches[0].replace('D','E').replace('d','e'))
        require(math.isfinite(row['gibbs_hartree']), 'Nonfinite Shermo Gibbs energy')
        row['shermo_log'] = str(log)
    low = min(r['gibbs_hartree'] for r in records)
    factors = [math.exp(-(r['gibbs_hartree']-low)*627.509474/(.00198720425864083*sher['T'])) for r in records]
    total = sum(factors)
    for row, factor in zip(records, factors):
        row['weight'] = factor/total
    return sorted(records, key=lambda r: r['gibbs_hartree'])


def run_candidate(candidate, cfg, run):
    from sampling import sample
    source = Path(candidate['input']).resolve()
    require(source.exists(), 'Missing input '+str(source))
    run.mkdir(parents=True, exist_ok=True)
    state = run/'state.json'
    identity = {'config': calculation_config(cfg), 'candidate': candidate, 'input_sha256': digest(source),
                'resources_sha256': {p.name:digest(p) for p in (ROOT/'resources').iterdir() if p.is_file()}}
    if state.exists():
        require(json.loads(state.read_text())['identity'] == identity, 'Run configuration changed; use new run_id')
    save(state, {'identity': identity, 'status': 'running'})
    ensemble = source if source.suffix.lower() == '.xyz' else sample(source, cfg, run/'sampling')
    blocks = xyz_read(ensemble)
    stereo = json.loads((run/'sampling/stereocenters.json').read_text()) if (run/'sampling/stereocenters.json').exists() else []
    limit = cfg['max_conformers']
    require(limit == 0 or limit > 0, 'Invalid conformer cap')
    selected = select_ensemble(blocks,limit)
    electrons = sum(ELEMENTS.index(a[0]) for a in selected[0][2])-cfg['charge']
    require(electrons >= cfg['multiplicity']-1 and (electrons-(cfg['multiplicity']-1)) % 2 == 0,
            'Electron count incompatible with charge/multiplicity')
    opt_rows, exclusions = [], []
    for index, comment, geom in selected:
        cid = 'c{:04d}'.format(index)
        work = run/'conformers'/cid
        work.mkdir(parents=True, exist_ok=True)
        source_name = 'Geom_source_{}_{}.xyz'.format(candidate['id'],cid.upper()) if cfg.get('output_layout')=='direct_v1' else 'source.xyz'
        (work/source_name).write_text(xyz_text(geom, comment))
        optimized = None
        for attempt in (1, 2, 3):
            attempt_dir = work/('opt'+str(attempt))
            try:
                optimized = quantum(ROOT/'resources'/('03_tem'+str(attempt)+'.gjf'), geom, cfg, attempt_dir, 'opt')
                if stereo:
                    from sampling import filter_stereochemistry
                    filter_stereochemistry([('',optimized['geometry'])],stereo)
                break
            except ValueError as exc:
                if (attempt_dir/'receipt.json').exists():
                    raise
                exclusions.append({'id': cid, 'attempt': attempt, 'reason': str(exc)})
                # A failed optimization remains archived. Retry starts from source geometry.
        require(optimized is not None, 'All optimization attempts failed for '+cid)
        opt_rows.append({'id': cid, 'geometry': optimized['geometry'], 'opt_energy': optimized['energy']})
    dedup = run/'dft_dedup'
    dedup.mkdir(exist_ok=True)
    opt_xyz = dedup/'optimized.xyz'
    opt_xyz.write_text(''.join(xyz_text(r['geometry'], str(r['opt_energy'])) for r in opt_rows))
    require(not (dedup/'cluster.xyz').exists() or (run/'conformers.json').exists(), 'Incomplete clustering stage')
    # isostat is rerun in a new private directory on resume; no stale generic output.
    cluster_dir = Path(tempfile.mkdtemp(prefix='cluster_',dir=str(dedup)))
    command([ROOT/'resources/isostat', opt_xyz, *cluster_options(cfg,'DFT'),
             '-nt', cfg['cores']], cluster_dir, cluster_dir/'isostat.log')
    kept, matched = [], set()
    for _, geom in xyz_read(cluster_dir/'cluster.xyz'):
        distances = sorted((rms(geom, r['geometry']), r['id'], r) for r in opt_rows)
        distance, cid, row = distances[0]
        require(distance < 1e-4 and cid not in matched, 'Cannot uniquely map DFT cluster to parent ID')
        matched.add(cid)
        kept.append(row)
    exclusions += [{'id': r['id'], 'reason': 'DFT isostat duplicate/threshold exclusion'} for r in opt_rows if r['id'] not in matched]
    records = []
    for row in kept:
        work = run/'conformers'/row['id']
        try:
            freq = quantum(ROOT/'resources/04_tem.gjf', row['geometry'], cfg, work/'freq', 'freq')
        except ValueError as exc:
            if 'imaginary/zero frequencies' not in str(exc):
                raise
            exclusions.append({'id': row['id'], 'reason': str(exc)})
            continue
        freq_dir = work/'freq'
        row.update(frequency=str(freq_dir/(artifact_stem(freq_dir,cfg,'Gau','freq')+'.log')), geometry=freq['geometry'])
        row['sp_energy'] = orca_sp(row['geometry'], cfg, work/'sp')
        records.append(row)
    require(bool(records), 'No frequency-qualified conformers')
    records = populations(records, cfg, run/'shermo')
    chosen, covered = [], 0.0
    for row in records:
        chosen.append(row)
        covered += row['weight']
        if covered >= cfg['population_coverage']:
            break
    for row in chosen:
        row['nmr_weight'] = row['weight']/covered
        row['nmr'] = {}
        for method in cfg['methods']:
            directory = run/'conformers'/row['id']/('nmr'+method)
            quantum(ROOT/'resources'/('06_tem'+template_method(method,cfg)+'_nmr.gjf'), row['geometry'], cfg, directory, 'nmr')
            row['nmr'][method] = str(directory/(artifact_stem(directory,cfg,'Gau','nmr')+'.log'))
    save(run/'conformers.json', {'candidate': candidate['id'], 'config': cfg, 'ensemble_total': len(blocks),
         'dft_selected': len(selected), 'coverage_within_computed_ensemble': covered,
         'records': records, 'nmr_ids': [r['id'] for r in chosen], 'exclusions': exclusions})
    save(state, {'identity': identity, 'status': 'complete'})


def validate_config(cfg):
    require(cfg.get('output_layout') in (None,'direct_v1'), 'Unknown output layout')
    require(all(isinstance(cfg[k],int) and not isinstance(cfg[k],bool)
                for k in ['charge','multiplicity','cores','memory_gb','orca_maxcore_mb','max_conformers']),
            'Charge, multiplicity, resources and conformer cap must be integers')
    require(cfg['cores'] > 0 and cfg['memory_gb'] > 0, 'Invalid resources')
    require(cfg['max_conformers'] >= 0 and cfg['orca_maxcore_mb'] > 0, 'Invalid conformer cap/MaxCore')
    require(isinstance(cfg.get('md_threads',4),int) and 0 < cfg.get('md_threads',4) <= cfg['cores'], 'Invalid MD thread count')
    require(cfg['orca_maxcore_mb']*cfg['cores'] < cfg['memory_gb']*1024*.85, 'ORCA MaxCore leaves insufficient memory margin')
    require(0 < cfg['population_coverage'] <= 1, 'Coverage must be (0,1]')
    require(cfg['methods'] and set(cfg['methods']) <= {'01','02','03'}, 'Unknown method')
    require(len(set(cfg['methods'])) == len(cfg['methods']), 'Duplicate NMR method')
    require(cfg['charge'] == int(cfg['charge']) and cfg['multiplicity'] >= 1, 'Invalid charge/spin')
    require(cfg['solvent'] == 'chloroform', 'Current calibrated interfaces support chloroform; version other solvents explicitly')
    require(cfg['shermo']['T'] > 0 and cfg['shermo']['P'] > 0,'Invalid thermochemistry state')
    require(cfg['candidates'], 'No candidates supplied')
    require(len(cfg['md_seeds']) == 3 and len(set(cfg['md_seeds'])) == 3
            and all(isinstance(v,int) and not isinstance(v,bool) and -2147483648 <= v <= 2147483647 and v != -1
                    for v in cfg['md_seeds']), 'Three distinct fixed MD seeds required')
    require(isinstance(cfg['grompp_maxwarn'],int) and cfg['grompp_maxwarn'] >= 0, 'Invalid grompp warning limit')
    require(0 <= cfg['resp2']['delta'] <= 1 and cfg['resp2']['solvent'] == 'water', 'RESP2 requires delta in [0,1] and water')
    validate_choices(cfg)
    parameters = json.loads((ROOT/'calibration.json').read_text())
    branches = set()
    for branch in cfg['postprocessing']:
        kind, method = branch['kind'], branch['method']
        require(kind in ('dp4','dp4plus','sts') and method in cfg['methods'], 'Unknown/disabled postprocessing branch')
        require((kind,method) not in branches, 'Duplicate postprocessing output')
        branches.add((kind,method))
        group = {'dp4':'dp4_classic','dp4plus':'dp4','sts':'sts'}[kind]
        require(branch['parameter_key'] in parameters[group], 'Unknown calibration parameter key')
        require(isinstance(branch.get('allow_exploratory',False),bool), 'allow_exploratory must be a boolean')
        if kind == 'dp4plus':
            expected = {mid:f+'|PCM|'+b for mid,(f,b,k) in method_map(cfg).items() if k=='dp4plus'}
            require(branch['parameter_key'] == expected.get(method), 'DP4+ parameter does not match this NMR template')
        if kind == 'sts':
            require(method_map(cfg)[method][2] == 'sts' and branch['parameter_key'] == 'CDCl3', 'STS requires the wB97XD branch and CDCl3 parameters')
    for c in cfg['candidates']:
        require(re.fullmatch(r'[A-Za-z0-9_-]+', c['id']) is not None, 'Unsafe candidate ID')
    require(len(set(c['id'] for c in cfg['candidates'])) == len(cfg['candidates']), 'Duplicate candidate ID')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['check','run','post'])
    parser.add_argument('config')
    args = parser.parse_args()
    config_path = Path(args.config).resolve()
    cfg = json.loads(config_path.read_text())
    for c in cfg['candidates']:
        c['input'] = str((config_path.parent/c['input']).resolve())
    cfg['output'] = str((config_path.parent/cfg['output']).resolve())
    if cfg.get('output_layout')=='direct_v1':
        require(cfg['output']==str(config_path.parent/'过程文件') and not (config_path.parent/'过程文件').is_symlink(),
                'Direct layout requires a real project/过程文件 directory')
    validate_config(cfg)
    if args.action == 'check':
        for name, path in cfg['software'].items():
            resolved = path if Path(path).is_absolute() else shutil.which(path)
            require(resolved and Path(resolved).exists(), 'Missing {}: {}'.format(name,path))
            if name != 'openmpi':
                require(os.access(resolved,os.X_OK), 'Software not executable: '+path)
        import numpy
        import openpyxl
        for filename in ['03_tem1.gjf','03_tem2.gjf','03_tem3.gjf','04_tem.gjf','05_tem.inp',
                         '06_tem01_nmr.gjf','06_tem02_nmr.gjf','06_tem03_nmr.gjf','isostat','sobtop','atomtype','shermo_settings.ini']:
            require((ROOT/'resources'/filename).exists(), 'Missing workflow resource '+filename)
        for c in cfg['candidates']:
            require(Path(c['input']).exists(), 'Missing input '+c['input'])
            if c['input'].lower().endswith('.xyz'):
                select_ensemble(xyz_read(c['input']),cfg['max_conformers'])
            else:
                for filename in ['md01.mdp','md02.mdp','md03.mdp','sobtop.ini','mol2_xyz','xyz2QC','bonded_param.dat',
                                 'assign_AT.dat','ATOMTYPE_GFF.DEF','LJ_param.dat']:
                    require((ROOT/'resources'/filename).exists(), 'Missing sampling resource '+filename)
                from sampling import md_parameters
                for i,seed in enumerate(cfg['md_seeds'],1):
                    md_parameters((ROOT/'resources'/('md%02d.mdp'%i)).read_text(),seed)
        print('Configuration and input checks passed; no computation submitted')
        return
    output = Path(cfg['output'])
    with output_lock(output):
        execute(args.action,cfg,output,config_path.parent)


def execute(action,cfg,output,config_dir):
    if action == 'run':
        require(bool(os.environ.get('PBS_JOBID')) or cfg.get('allow_local_compute') is True,
                'Computations require a PBS allocation')
        if os.environ.get('PBS_JOBID'):
            nodefile = Path(os.environ.get('PBS_NODEFILE','/nonexistent-pbs-nodefile'))
            require(nodefile.is_file(), 'PBS_NODEFILE is required to verify the actual single-node allocation')
            slots = [line.strip() for line in nodefile.read_text().splitlines() if line.strip()]
            require(slots and len(set(slots)) == 1 and cfg['cores'] <= len(slots),
                    'This workflow requires enough allocated CPU slots on one node')
            if Path('/proc/meminfo').exists():
                available = re.search(r'^MemAvailable:\s+(\d+)',Path('/proc/meminfo').read_text(),re.M)
                require(available and int(available[1]) > (cfg['memory_gb']+4)*1024*1024,
                        'Node available memory is below application memory plus safety margin')
        for candidate in cfg['candidates']:
            run_candidate(candidate, cfg, output/candidate['id'])
    from postprocess import export
    export(cfg, output, config_dir)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, subprocess.SubprocessError) as exc:
        sys.exit('STOP: '+str(exc))
