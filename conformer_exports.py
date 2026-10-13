"""Coordinates and auditable population summaries, independent of NMR properties."""
import csv
import json
import math
from pathlib import Path
import tempfile
import time
from module_options import analysis


def electronic_records(blocks,temperature):
    from nmrflow import ensemble_energy,require
    require(temperature>0,'Invalid population temperature')
    rows=[{'id':'c{:04d}'.format(i),'geometry':g,'electronic_energy_hartree':ensemble_energy(c)} for i,(c,g) in enumerate(blocks,1)]
    low=min(r['electronic_energy_hartree'] for r in rows)
    factors=[math.exp(-(r['electronic_energy_hartree']-low)*627.509474/(.00198720425864083*temperature)) for r in rows]
    for r,v in zip(rows,factors):r['weight']=v/sum(factors)
    return sorted(rows,key=lambda r:(-r['weight'],r['id']))


def export_candidate(manifest,destination,top_n):
    from nmrflow import require,xyz_text,save
    records=manifest['records']
    require(records and all(math.isfinite(r['weight']) and r['weight']>=0 for r in records),'Invalid conformer weights')
    require(abs(sum(r['weight'] for r in records)-1)<1e-8,'Conformer populations do not sum to one')
    records=sorted(records,key=lambda r:(-r['weight'],r['id']))
    selected=records[:top_n];coverage=sum(r['weight'] for r in selected)
    key='electronic_energy_hartree' if manifest.get('weight_model')=='gfn2_electronic_energy' else 'gibbs_hartree'
    destination=Path(destination);destination.mkdir(parents=True,exist_ok=True)
    lowest=min(r[key] for r in records);ids={r['id'] for r in selected}
    headers=['conformer','weight_model','relative_energy_kcal_mol','full_ensemble_percent','exported','selected_normalized_percent','selected_coverage_percent']
    with (destination/'populations.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=headers);w.writeheader()
        for r in records:w.writerow(dict(zip(headers,[r['id'],key,(r[key]-lowest)*627.509474,100*r['weight'],r['id'] in ids,100*r['weight']/coverage if r['id'] in ids else '',100*coverage])))
    texts=[]
    for r in selected:
        text=xyz_text(r['geometry'],'{}; population={:.12f}; selected_share={:.12f}; model={}'.format(r['id'],r['weight'],r['weight']/coverage,key))
        (destination/(r['id']+'.xyz')).write_text(text,encoding='utf-8');texts.append(text)
    (destination/'top_conformers.xyz').write_text(''.join(texts),encoding='utf-8')
    summary={'candidate':manifest['candidate'],'weight_model':key,'requested_top_n':top_n,'exported':len(selected),'ensemble_size':len(records),'selected_coverage':coverage,
             'ids':[r['id'] for r in selected],'temperature_K':manifest['config']['shermo']['T'],
             'scope':'GFN2 sampled cluster only' if key=='electronic_energy_hartree' else 'frequency-qualified computed DFT subset only',
             'warning':'Electronic-energy estimates omit thermal/entropy corrections' if key=='electronic_energy_hartree' else 'Finite sampling and DFT budget; not total conformational-space coverage'}
    save(destination/'summary.json',summary)
    return summary


def export_all(cfg,output):
    from nmrflow import require,calculation_config,save,digest
    parent=Path(output)/'module_exports';parent.mkdir(parents=True,exist_ok=True)
    pending=Path(tempfile.mkdtemp(prefix='.conformers_pending_',dir=str(parent)))
    result={}
    for c in cfg['candidates']:
        path=Path(output)/c['id']/'conformers.json'
        m=json.loads(path.read_text(encoding='utf-8'))
        require(calculation_config(m['config'])==calculation_config(cfg),'Conformer manifest/config mismatch')
        result[c['id']]=export_candidate(m,pending/c['id'],analysis(cfg)['top_n'])
        result[c['id']]['manifest_sha256']=digest(path)
    save(pending/'summary.json',result)
    final=parent/'conformers'
    if final.exists():final.rename(parent/('conformers_previous_'+time.strftime('%Y%m%d_%H%M%S')+'_'+pending.name[-6:]))
    pending.rename(final)
