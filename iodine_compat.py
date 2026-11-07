"""Explicit iodine basis compatibility; not a claim of spin-orbit NMR accuracy."""
from pathlib import Path
import re

POLICY='iodine_mixed_v1'
NOTE=('Iodine compatibility: light atoms retain the original basis; I uses def2-SVP/ECP for RESP/geometry/frequency, '
      'all-electron DGDZVP for legacy Gaussian NMR. No spin-orbit shielding is included; C-I signals require separate review. '
      'Original DP4+/STS calibrations are not validated for this mixed-basis iodine protocol.')


def input_elements(path):
    path=Path(path);text=path.read_text(encoding='utf-8-sig')
    if path.suffix.lower()=='.pdb':
        result=[line[76:78].strip().title() for line in text.splitlines() if line.startswith(('ATOM  ','HETATM'))]
    else:
        result=[]
        for line in text.splitlines():
            fields=line.split()
            if len(fields)!=4 or not re.fullmatch(r'[A-Z][a-z]?',fields[0]):continue
            try:[float(v) for v in fields[1:]]
            except ValueError:continue
            result.append(fields[0])
    if not result or any(not x for x in result):raise ValueError('Cannot determine explicit input elements: '+str(path))
    return result


def configure(cfg,elements):
    if 'I' not in elements:return cfg
    if cfg.get('element_policy') not in (None,POLICY):raise ValueError('Unknown iodine element policy')
    cfg['element_policy']=POLICY
    cfg['element_policy_elements']=sorted(set(elements))
    cfg['iodine_nmr_scope']='basis-compatible scalar/nonrelativistic legacy preview; spin-orbit correction absent'
    for branch in cfg.get('postprocessing',[]):
        if NOTE not in branch['differences']:branch['differences'].append(NOTE)
        branch['allow_exploratory']=False
    return cfg


def configure_inputs(cfg,paths):
    elements=[]
    for path in paths:elements.extend(input_elements(path))
    return configure(cfg,elements)


def require_policy(cfg,elements):
    if 'I' in elements and cfg.get('element_policy')!=POLICY:
        raise ValueError('Input contains iodine, but the default 6-31G/6-311G templates do not cover I. '
                         'Create a new iodine-compatible project; do not submit the unmodified legacy templates.')


def render(text,elements,cfg,stage):
    if 'I' not in elements:return text
    require_policy(cfg,elements)
    # Explicit custom basis decks must not be silently replaced.
    routes=re.findall(r'(?m)^#[^\n]*$',text)
    if len(routes)!=1:raise ValueError('Iodine adapter requires one explicit Gaussian route')
    route=routes[0]
    match=re.search(r'/((?:6-31(?:1)?)(?:\+\+|\+)?G(?:\*\*|\*|\([^\s]+\))?)',route,re.I)
    if match is None:raise ValueError('Iodine adapter only replaces verified legacy Pople templates; custom inputs require review')
    non_i=sorted(set(elements)-{'I'})
    if any(x not in {'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr'} for x in non_i):
        raise ValueError('Additional unsupported heavy element; iodine-only adapter is insufficient')
    nmr=stage in ('nmr','coupling','coupling_shielding','ecd')
    use_ecp=not nmr
    changed=route[:match.start()]+'/'+('GenECP' if use_ecp else 'Gen')+route[match.end():]
    if stage=='resp':
        changed,count=re.subn(r'(?i)\bpop\s*=\s*MK\b','Pop=(MK,ReadRadii)',changed)
        if count!=1:raise ValueError('Iodine RESP adapter requires the verified Pop=MK template')
    if not re.search(r'(?i)\b(?:5d|6d)\b',changed):changed+=' 6D'
    text=text.replace(route,changed,1)
    tail=[]
    if non_i:tail+=[' '.join(non_i)+' 0',match[1],'****']
    tail+=['I 0','Def2SVP' if use_ecp else 'DGDZVP','****','']
    if use_ecp:tail+=['I 0','Def2SVP','']
    # Bondi iodine radius in Angstrom; Gaussian has no default MK radius for I.
    if stage=='resp':tail+=['I 1.98','']
    return text.rstrip()+'\n\n'+'\n'.join(tail)+'\n'
