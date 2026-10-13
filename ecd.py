"""TDDFT ECD transitions, Boltzmann averaging and traceable unaltered plots."""
import csv
import json
import math
from pathlib import Path
import re
import tempfile
import time
from module_options import analysis,options

HC_EV_NM=1239.8419843320026


def parse_transitions(text):
    from nmrflow import require,NUMBER
    states={}
    pattern=r'Excited State\s+(\d+):[^\n]*?('+NUMBER+r')\s+eV\s+('+NUMBER+r')\s+nm\s+f=\s*('+NUMBER+r')'
    for n,e,nm,f in re.findall(pattern,text):
        n=int(n)
        require(n not in states,'Multiple TD excited-state tables; split jobs first')
        states[n]={'state':n,'energy_ev':float(e),'wavelength_nm':float(nm),'oscillator_strength':float(f)}
    require(states and set(states)==set(range(1,len(states)+1)),'Missing/noncontiguous TD states')
    for gauge in ('length','velocity'):
        headers=list(re.finditer(r'^\s*state[^\n]*R\('+gauge+r'\)[^\n]*\n',text,re.M))
        require(len(headers)==1,'Expected one ECD R('+gauge+') table')
        seen=set()
        for line in text[headers[0].end():].splitlines():
            p=line.split()
            if len(p)<5 or not re.fullmatch(r'\d+',p[0]):break
            n=int(p[0]);require(n in states and n not in seen,'Invalid ECD rotatory-strength state')
            value=float(p[4].replace('D','E').replace('d','e'))
            require(math.isfinite(value),'Nonfinite ECD rotatory strength')
            states[n]['r_'+gauge]=value;seen.add(n)
        require(seen==set(states),'Incomplete ECD rotatory-strength table')
    for s in states.values():
        require(all(math.isfinite(v) for v in s.values()) and s['energy_ev']>0 and s['wavelength_nm']>0,'Invalid TD transition')
        require(abs(s['energy_ev']*s['wavelength_nm']-HC_EV_NM)<.2,'TD energy/wavelength mismatch')
    return [states[n] for n in sorted(states)]


def spectrum(transitions,wavelengths,fwhm_ev,gauge='velocity',energy_scale=1.0,shift_nm=0.0):
    """Gaussian broadening, Δε in L mol^-1 cm^-1; Gaussian R unit is 10^-40 cgs.

    FWHM=2 sqrt(ln 2) width, where width is the half width at 1/e height.
    Nugroho & Morita, DOI 10.1007/s11418-013-0768-x, equation 1.
    Energy and Gaussian width are both in eV; their conversion factors cancel.
    A wavelength shift is a display translation, not a fitted correction.
    """
    import numpy as np
    from nmrflow import require
    require(fwhm_ev>0 and energy_scale>0 and gauge in ('length','velocity'),'Invalid ECD broadening')
    nm=np.asarray(wavelengths,dtype=float)-shift_nm
    require(np.all(nm>0),'ECD shift produces nonpositive wavelength')
    width=fwhm_ev/(2*math.sqrt(math.log(2)))
    energy=HC_EV_NM/nm;curve=np.zeros_like(energy)
    for t in transitions:
        ei=t['energy_ev']*energy_scale
        curve+=ei*t['r_'+gauge]*np.exp(-((energy-ei)/width)**2)
    return curve/(22.96*width*math.sqrt(math.pi))


def experimental(path):
    from nmrflow import require
    from postprocess import table
    if not Path(path).is_file():return None
    rows=table(path);require(len(rows)>=2,'ECD experimental curve requires at least two points')
    values=[(float(r['wavelength_nm']),float(r['delta_epsilon'])) for r in rows]
    require(all(math.isfinite(x) and math.isfinite(y) and x>0 for x,y in values),'Invalid experimental ECD point')
    require(len({x for x,y in values})==len(values),'Duplicate experimental wavelength')
    return values  # Do not truncate, smooth, clip, replace a baseline or change signs.


def plot(path,wavelengths,curves,exp=None,normalize=False,title='ECD',individual=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    fig,ax=plt.subplots(figsize=(8,4.8),constrained_layout=True)
    def scaled(y):
        y=np.asarray(y,dtype=float);peak=np.max(np.abs(y))
        return y/peak if normalize and peak else y
    if individual and not normalize:
        for label,y in individual.items():ax.plot(wavelengths,y,color='#9aa7b0',alpha=.35,lw=.65)
    for label,y in curves.items():ax.plot(wavelengths,scaled(y),lw=1.65,label=label)
    if exp is not None:ax.plot([r[0] for r in exp],scaled([r[1] for r in exp]),color='black',lw=1.8,label='Experiment (raw)')
    ax.axhline(0,color='#77848d',lw=.7);ax.set(xlabel='Wavelength / nm',ylabel='Normalized ECD / a.u.' if normalize else r'$\Delta\epsilon$ / L mol$^{-1}$ cm$^{-1}$',title=title)
    ax.legend(frameon=False);ax.grid(alpha=.15)
    for suffix in ('.png','.svg'):fig.savefig(Path(path).with_suffix(suffix),dpi=200)
    plt.close(fig)


def export_all(cfg,output,config_dir):
    import numpy as np
    from nmrflow import require,save,digest,gaussian_read,rms,calculation_config
    from postprocess import write_csv
    output,config_dir=Path(output),Path(config_dir);a,m=analysis(cfg),options(cfg)
    wavelengths=np.arange(a['ecd_min_nm'],a['ecd_max_nm']+a['ecd_step_nm']*.5,a['ecd_step_nm'])
    exp_path=config_dir/a['ecd_experiment'];exp=experimental(exp_path)
    parent=output/'module_exports';parent.mkdir(exist_ok=True)
    pending=Path(tempfile.mkdtemp(prefix='.ecd_pending_',dir=str(parent)))
    allcurves={mid:{} for mid in m['ecd_methods']};provenance=[];coverage={}
    for candidate in cfg['candidates']:
        cid=candidate['id'];manifest=output/cid/'conformers.json';state=json.loads(manifest.read_text())
        require(calculation_config(state['config'])==calculation_config(cfg),'ECD manifest/config mismatch')
        chosen=[r for r in state['records'] if r['id'] in state['nmr_ids']]
        require(chosen and abs(sum(r['nmr_weight'] for r in chosen)-1)<1e-8,'Invalid ECD weights')
        coverage[cid]=state['coverage_within_computed_ensemble']
        for mid in m['ecd_methods']:
            dest=pending/cid/('method'+mid);dest.mkdir(parents=True)
            raw_curve=np.zeros_like(wavelengths);display=np.zeros_like(wavelengths);individual={};transitions=[]
            for row in chosen:
                log=Path(row['ecd'][mid]);receipt=json.loads((log.parent/'receipt.json').read_text())
                require(receipt['output_sha256']==digest(log) and receipt['input_sha256']==digest(log.with_suffix('.gjf')) and receipt['stage']=='ecd','ECD receipt/input/output hash mismatch')
                parsed=gaussian_read(log,'ecd')
                require(rms(parsed['geometry'],row['geometry'])<1e-4,'ECD geometry/manifest mismatch')
                require(len(parsed['transitions'])==m['ecd_states'],'Unexpected ECD state count')
                provenance.append({'candidate':cid,'method':mid,'conformer':row['id'],'sha256':digest(log),'path':str(log),'weight':row['nmr_weight']})
                for t in parsed['transitions']:transitions.append(dict(t,conformer=row['id'],weight=row['nmr_weight']))
                curve=spectrum(parsed['transitions'],wavelengths,a['ecd_fwhm_ev'],a['ecd_gauge'])
                raw_curve+=row['nmr_weight']*curve;individual[row['id']]=curve
                display+=row['nmr_weight']*spectrum(parsed['transitions'],wavelengths,a['ecd_fwhm_ev'],a['ecd_gauge'],a['ecd_energy_scale'],a['ecd_shift_nm'])
            write_csv(dest/'transitions.csv',list(transitions[0]),transitions)
            headers=['wavelength_nm','raw_delta_epsilon','display_delta_epsilon']+list(individual)
            with (dest/'curves.csv').open('w',encoding='utf-8-sig',newline='') as f:
                w=csv.writer(f);w.writerow(headers)
                w.writerows([[nm,raw_curve[i],display[i]]+[v[i] for v in individual.values()] for i,nm in enumerate(wavelengths)])
            plot(dest/'ecd_raw',wavelengths,{'Boltzmann average (raw)':raw_curve},exp,title=cid+' / ECD '+mid,individual=individual)
            adjusted=a['ecd_energy_scale']!=1.0 or a['ecd_shift_nm']!=0
            if adjusted:plot(dest/'ecd_display',wavelengths,{'Raw':raw_curve,'Explicitly adjusted':display},exp,title=cid+' / explicit display correction')
            if a['ecd_normalize']:plot(dest/'ecd_normalized',wavelengths,{'Boltzmann average':display},exp,normalize=True,title=cid+' / independently normalized curves')
            allcurves[mid][cid]=raw_curve
    for mid,curves in allcurves.items():plot(pending/('all_candidates_method'+mid),wavelengths,curves,exp,title='ECD method '+mid+' / raw curves')
    if exp is not None:
        import shutil
        shutil.copy2(exp_path,pending/'experiment_original.csv')
    save(pending/'metadata.json',{'status':'ecd_exports_complete','settings':a,'modules':m,'coverage':coverage,'provenance':provenance,
         'experimental_sha256':digest(exp_path) if exp is not None else None,
         'units':'delta_epsilon: L mol^-1 cm^-1; Gaussian R: 10^-40 cgs',
         'protocol_difference':'Gas-phase optimized geometries and gas-phase ORCA/Shermo populations; SMD solvent used in TDDFT only',
         'broadening_reference':'https://doi.org/10.1007/s11418-013-0768-x',
         'notes':['No experimental clipping/smoothing/baseline replacement/mirroring','No automatic shift, scale or candidate probability fit','30 states may not cover the requested short-wavelength range; inspect transitions and convergence']})
    final=parent/'ecd'
    if final.exists():final.rename(parent/('ecd_previous_'+time.strftime('%Y%m%d_%H%M%S')+'_'+pending.name[-6:]))
    pending.rename(final)
