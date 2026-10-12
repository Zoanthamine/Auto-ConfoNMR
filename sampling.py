"""Original MD -> VMD -> GFN1/GFN2 sampling, isolated and checked."""
import json
from pathlib import Path
import re
import shutil
import math
import os

from nmrflow import ROOT, command, require, save, xyz_read, xyz_text, geometry_text, digest


def stereo_volume(geom, neighbors):
    a,b,c,d = [geom[i-1][1:] for i in neighbors]
    x,y,z = [[p[j]-d[j] for j in range(3)] for p in (a,b,c)]
    return x[0]*(y[1]*z[2]-y[2]*z[1])-x[1]*(y[0]*z[2]-y[2]*z[0])+x[2]*(y[0]*z[1]-y[1]*z[0])


def pdb_stereocenters(source, geom):
    lines = source.read_text().splitlines()
    atoms = [s for s in lines if s.startswith(('ATOM  ','HETATM'))]
    serial_map = {int(s[6:11]):i for i,s in enumerate(atoms,1)}
    bonds = {i:set() for i in range(1,len(geom)+1)}
    for line in lines:
        if line.startswith('CONECT'):
            ids = [int(line[i:i+5]) for i in range(6,len(line),5) if line[i:i+5].strip()]
            for value in ids[1:]:
                a,b = serial_map[ids[0]],serial_map[value]
                bonds[a].add(b);bonds[b].add(a)
    reference = []
    for i,neighbors in bonds.items():
        # Conservative guard for tetrahedral carbon centers; identical-H centers
        # are excluded. This is not a full CIP assignment or graph-isomorphism test.
        if geom[i-1][0] == 'C' and len(neighbors) == 4 and sum(geom[n-1][0]=='H' for n in neighbors) <= 1:
            colors = {j:geom[j-1][0] for j in bonds if j != i}
            for _ in range(len(geom)):
                keys = {j:(colors[j],tuple(sorted(colors[k] for k in bonds[j] if k != i))) for j in colors}
                ranks = {k:str(n) for n,k in enumerate(sorted(set(keys.values())))}
                colors = {j:ranks[k] for j,k in keys.items()}
            if len({colors[j] for j in neighbors}) != 4:
                continue
            ns = sorted(neighbors)
            value = stereo_volume(geom,ns)
            require(abs(value) > .1, 'Nearly planar starting stereocenter '+str(i))
            reference.append({'center':i,'neighbors':ns,'sign':1 if value>0 else -1})
    return reference


def filter_stereochemistry(blocks, reference):
    kept,rejected = [],[]
    for index,(comment,geom) in enumerate(blocks,1):
        bad = [r['center'] for r in reference if stereo_volume(geom,r['neighbors'])*r['sign'] <= .1]
        if bad: rejected.append({'block':index,'centers':bad})
        else: kept.append((comment,geom))
    require(kept, 'No conformers preserve the input stereochemistry')
    return kept,rejected


def md_parameters(text, seed):
    """Require a full annealing cycle and fix both velocity and thermostat seeds."""
    text, count = re.subn(r'^[ \t]*gen_seed[ \t]*=.*$', 'gen_seed = '+str(seed), text, flags=re.M)
    require(count == 1, 'MD seed setting missing/duplicated')
    text = re.sub(r'(?m)^[ \t]*ld[-_]seed[ \t]*=.*\n?', '', text)
    text += '\nld_seed = '+str(seed % 2147483647)+'\n'
    def setting(name):
        values = re.findall(r'(?mi)^[ \t]*'+name+r'[ \t]*=[ \t]*([^;\n]+)',text)
        require(len(values) == 1, 'Missing/duplicate MD parameter '+name)
        return values[0].strip()
    steps,dt = int(setting('nsteps')),float(setting('dt'))
    times = [float(v) for v in setting('annealing_time').split()]
    temps = [float(v) for v in setting('annealing_temp').split()]
    require(steps > 0 and dt > 0 and len(times) == len(temps) == int(setting('annealing_npoints')), 'Invalid annealing schedule')
    require(times[0] == 0 and all(a < b for a,b in zip(times,times[1:])), 'Invalid annealing time order')
    require(steps*dt >= times[-1] and abs(temps[0]-temps[-1]) < 1e-8, 'MD must complete a full cooling cycle')
    return text, steps*dt


def resp2(mol2, cfg, work):
    """Two real ESP/RESP fits; delta weights solvent phase, not absent columns."""
    software = cfg['software']
    command([software['obabel'], mol2, '-O', work/'input.xyz'], work, work/'conversion.log')
    blocks = xyz_read(work/'input.xyz')
    require(len(blocks) == 1, 'Sampling input must contain exactly one molecule')
    geom = blocks[0][1]
    delta = cfg['resp2']['delta']
    require(0 <= delta <= 1, 'Invalid RESP2 delta')
    values = {}
    for phase in ['gas','solvent']:
        directory = work/('resp_'+phase)
        directory.mkdir(exist_ok=True)
        route = '# B3LYP/6-31g pop=MK IOp(6/33=2,6/42=6)'
        if phase == 'solvent':
            route += ' scrf(solvent={})'.format(cfg['resp2']['solvent'])
        stem = 'Gau_RESP_{}_{}'.format(phase,work.parent.name) if cfg.get('output_layout')=='direct_v1' else 'gau'
        inp = ('%mem={}GB\n%nprocshared={}\n%chk='+stem+'.chk\n{}\n\nRESP2 {}\n\n{} {}\n{}\n\n').format(
            cfg['memory_gb'], cfg['cores'], route, phase, cfg['charge'], cfg['multiplicity'], geometry_text(geom))
        (directory/(stem+'.gjf')).write_text(inp)
        env = os.environ.copy()
        env['NMR_DETAIL_LOG'] = '0'
        (directory/'scratch').mkdir(exist_ok=True)
        env['GAUSS_SCRDIR'] = str((directory/'scratch').resolve())
        command([software['g16']], directory, directory/(stem+'.out'), env=env, stdin=inp)
        out = (directory/(stem+'.out')).read_text(errors='replace')
        require('Normal termination of Gaussian' in out and 'Error termination' not in out, 'RESP2 Gaussian failed')
        command([software['formchk'], stem+'.chk', stem+'.fchk'], directory, directory/(stem+'_formchk.log'))
        controls = '7\n18\n8\n1\n'+stem+'.out\ny\n0\n0\nq\n'
        (directory/'multiwfn.in').write_text(controls)
        command([software['multiwfn'], stem+'.fchk', '-ispecial', '1'], directory, directory/(stem+'_Multiwfn.log'), stdin=controls)
        chg = [line.split() for line in (directory/(stem+'.chg')).read_text().splitlines() if line.strip()]
        require(len(chg) == len(geom) and all(len(x) == 5 for x in chg), 'RESP charge table incomplete')
        require([x[0] for x in chg] == [x[0] for x in geom], 'RESP atom order mismatch')
        values[phase] = [float(x[4]) for x in chg]
        require(all(math.isfinite(v) for v in values[phase]), 'Nonfinite RESP charge')
        require(abs(sum(values[phase])-cfg['charge']) < .002, 'RESP charge sum incorrect')
    charges = [(1-delta)*a+delta*b for a,b in zip(values['gas'], values['solvent'])]
    require(abs(sum(charges)-cfg['charge']) < .002, 'RESP2 total charge incorrect')
    (work/'molecule.chg').write_text('\n'.join('{} {:.8f} {:.8f} {:.8f} {:.10f}'.format(*(a+[q]))
                                           for a,q in zip(geom, charges))+'\n')
    save(work/'resp2.json', {'delta': delta, 'solvent': cfg['resp2']['solvent'],
         'protocol':'legacy-level two-phase RESP mixing; not the published RESP2 electronic-structure protocol',
         'level': 'B3LYP/6-31g', 'gas': values['gas'], 'solvent_charges': values['solvent'], 'mixed': charges})
    return geom


def sample(source, cfg, work):
    work.mkdir(parents=True, exist_ok=True)
    receipt, final = work/'receipt.json', work/'gfn2.xyz'
    if receipt.exists():
        saved = json.loads(receipt.read_text())
        require(saved['source_sha256'] == digest(source) and saved['output_sha256'] == digest(final), 'Sampling files changed')
        xyz_read(final)
        return final
    require(not (work/'molecule.mol2').exists(), 'Incomplete sampling; preserve old directory and use a new run_id')
    resources = ROOT/'resources'
    for p in resources.iterdir():
        if p.is_file():
            shutil.copy2(p, work/p.name)
    sw = cfg['software']
    command([sw['obabel'], source, '-O', 'molecule.mol2', '--partialcharge', 'none'], work, work/'obabel.log')
    geom = resp2(work/'molecule.mol2', cfg, work)
    stereo = []
    if source.suffix.lower() == '.pdb':
        atoms = [line for line in source.read_text().splitlines() if line.startswith(('ATOM  ','HETATM'))]
        require(len(atoms) == len(geom), 'PDB conversion changed atom count')
        original = [[line[76:78].strip()]+[float(line[a:b]) for a,b in [(30,38),(38,46),(46,54)]] for line in atoms]
        from nmrflow import rms
        require(rms(original,geom) < 1e-3, 'PDB conversion changed atom order/geometry')
        save(work/'atom_map.json', [{'input_serial':int(line[6:11]),'calculation_atom':i,'element':a[0]}
                                  for i,(line,a) in enumerate(zip(atoms,geom),1)])
        stereo = pdb_stereocenters(source,geom)
        save(work/'stereocenters.json',stereo)
    # Accept both .top and .itp default paths before returning to the main menu.
    controls = '7\n10\n\nmolecule.chg\n0\n2\n\n1\n2\n4\n\n\n0\n'
    (work/'sobtop.in').write_text(controls)
    command([work/'sobtop', 'molecule.mol2'], work, work/'sobtop.log', stdin=controls)
    require(all((work/name).is_file() for name in ['molecule.top','molecule.itp','molecule.gro']),
            'Sobtop missing topology/include/coordinates')
    topology_log = (work/'sobtop.log').read_text(errors='replace')
    topology_warnings = re.findall(r'Missing (?:dihedral|improper) types[^\n]*(?:\n[ \t]+\d+[^\n]*)*',topology_log)
    save(work/'topology_quality.json', {'protocol':'legacy prebuilt parameters with guessed missing terms',
         'warnings':topology_warnings,'fully_parameterized':not topology_warnings})
    command([sw.get('gmx','gmx'),'editconf','-f','molecule.gro','-o','boxed.gro','-d','2'], work, work/'editconf.log')
    complete, failures = [], []
    for index, seed in enumerate(cfg['md_seeds'], 1):
        prefix = 'md{:02d}'.format(index)
        text = (resources/(prefix+'.mdp')).read_text(errors='replace')
        text, duration = md_parameters(text,seed)
        (work/(prefix+'.mdp')).write_text(text)
        try:
            command([sw.get('gmx','gmx'),'grompp','-f',prefix+'.mdp','-c','boxed.gro','-p','molecule.top','-o',prefix+'.tpr',
                     '-maxwarn',cfg['grompp_maxwarn']],work,work/(prefix+'_grompp.log'))
            env = os.environ.copy()
            threads = min(cfg['cores'],cfg.get('md_threads',4))
            env['OMP_NUM_THREADS'] = str(threads)
            command([sw.get('gmx','gmx'),'mdrun','-deffnm',prefix,'-ntmpi','1','-ntomp',threads,'-nt',threads],
                    work,work/(prefix+'_run.log'),env=env)
            log = (work/(prefix+'.log')).read_text(errors='replace')+(work/(prefix+'_run.log')).read_text(errors='replace')
            require('Finished mdrun' in log and not re.search(r'LINCS WARNING|Fatal error|Segmentation fault|\bnan\b',log,re.I), 'MD incomplete/unstable')
            command([sw.get('gmx','gmx'),'check','-f',prefix+'.xtc'],work,work/(prefix+'_check.log'))
            check = (work/(prefix+'_check.log')).read_text(errors='replace')
            last = re.findall(r'Last frame\s+\d+\s+time\s+([\d.]+)',check)
            require(last and float(last[-1]) >= duration-.01, 'MD trajectory did not reach requested duration')
            complete.append(prefix)
        except (ValueError, OSError) as exc:
            failures.append({'trajectory':prefix,'reason':str(exc)})
    require(bool(complete), 'No complete stable MD trajectory')
    preferred = complete[-1]
    command([sw.get('gmx','gmx'),'trjconv','-s',preferred+'.tpr','-f',preferred+'.xtc','-o','whole.xtc','-pbc','whole'],
            work,work/'trjconv.log',stdin='0\n')
    # Three independent replicas; the verified route selects md03 -> md02 -> md01.
    tcl = ('mol new '+preferred+'.gro type gro waitfor all\n'
           'mol addfile whole.xtc type xtc waitfor all\n'
           'set sel [atomselect top all]\n'
           'animate write xyz trajectory.xyz beg 1 end -1 waitfor all sel $sel\nquit\n')
    (work/'extract.tcl').write_text(tcl)
    (work/'empty.vmdrc').write_text('')
    env = os.environ.copy()
    env['VMDDIR'] = str(Path(sw['vmd']).parent)
    command([sw['vmd'],'-dispdev','text','-startup','empty.vmdrc','-e','extract.tcl'],work,work/'vmd.log',env=env)
    trajectory = xyz_read(work/'trajectory.xyz')
    require([a[0] for a in trajectory[0][1]] == [a[0] for a in geom], 'VMD changed atom order')
    from protocol_options import cluster_options
    previous = work/'trajectory.xyz'
    stereo_rejected = {}
    if stereo:
        trajectory,rejected = filter_stereochemistry(trajectory,stereo)
        stereo_rejected['MD'] = rejected
        previous = work/'trajectory_stereo_checked.xyz'
        previous.write_text(''.join(xyz_text(g,c) for c,g in trajectory))
    for gfn in (1,2):
        stage = work/('gfn'+str(gfn))
        stage.mkdir()
        command([sw['crest'],'-mdopt',previous,'-chrg',cfg['charge'],'-uhf',cfg['multiplicity']-1,
                 '-gfn',gfn,'--screen','-opt','normal','-nproc',cfg['cores']],stage,stage/'crest.log')
        log = (stage/'crest.log').read_text(errors='replace')
        require(re.search(r'CREST terminated normally|normal termination of CREST',log,re.I) is not None,
                'CREST normal termination marker missing')
        conformers = xyz_read(stage/'crest_ensemble.xyz')
        require([a[0] for a in conformers[0][1]] == [a[0] for a in geom], 'CREST changed atom identities/order')
        if stereo:
            conformers,rejected = filter_stereochemistry(conformers,stereo)
            stereo_rejected['GFN'+str(gfn)] = rejected
        (stage/'stereo_checked.xyz').write_text(''.join(xyz_text(g,c) for c,g in conformers))
        command([resources/'isostat','stereo_checked.xyz',*cluster_options(cfg,'GFN'+str(gfn)),
                 '-nt',cfg['cores']],stage,stage/'isostat.log')
        xyz_read(stage/'cluster.xyz')
        previous = stage/'cluster.xyz'
    shutil.copy2(previous,final)
    save(work/'sampling.json', {'selected_trajectory':preferred,'complete':complete,'failures':failures,
         'topology_warnings':topology_warnings,
         'stereochemistry_rejections':stereo_rejected,
         'md_seeds':cfg['md_seeds'],'frames':len(trajectory),'atom_reference':'input.xyz',
         'note':'Input atom map refers to converted MOL2/XYZ order; verify it before experimental comparison.'})
    save(receipt, {'source_sha256':digest(source),'output_sha256':digest(final)})
    return final
