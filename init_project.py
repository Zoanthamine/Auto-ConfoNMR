#!/usr/bin/env python3
"""Create a project from explicitly supplied candidate structures."""
import argparse
import csv
import json
from pathlib import Path
import re
import shutil

ROOT = Path(__file__).resolve().parent


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory')
    p.add_argument('structures',nargs='+')
    p.add_argument('--charge',type=int,default=0)
    p.add_argument('--multiplicity',type=int,default=1)
    args = p.parse_args()
    target = Path(args.directory).resolve()
    if target.exists():
        p.error('Project directory already exists; choose a new directory')
    structures = [Path(x).resolve() for x in args.structures]
    if not all(x.is_file() for x in structures):
        p.error('One or more structure files do not exist')
    ids = [re.sub(r'[^A-Za-z0-9_-]','_',x.stem) for x in structures]
    if len(set(ids)) != len(ids):
        p.error('Candidate names are not unique after normalization')
    (target/'inputs').mkdir(parents=True)
    cfg = json.loads((ROOT/'job.example.json').read_text())
    cfg.update(charge=args.charge,multiplicity=args.multiplicity,output='runs/run_001',candidates=[])
    for cid,source in zip(ids,structures):
        name = cid+source.suffix.lower()
        shutil.copy2(source,target/'inputs'/name)
        cfg['candidates'].append({'id':cid,'input':'inputs/'+name})
    (target/'job.json').write_text(json.dumps(cfg,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    for name,columns in [('experiment',['signal','nucleus','ppm']),
                         ('assignments',['candidate','signal','atoms','hybridization','carbon_type','heavy_atom'])]:
        # Templates remain separate until real assignments are supplied. Empty tables
        # must not be mistaken for a completed experimental comparison.
        with (target/(name+'.template.csv')).open('w',encoding='utf-8-sig',newline='') as f:
            csv.writer(f).writerow(columns)
    print('Project prepared: '+str(target/'job.json'))
    print('Use submit.sh with that configuration. Fill/rename experimental templates when assignments are available.')


if __name__ == '__main__':
    main()
