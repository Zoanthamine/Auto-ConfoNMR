"""Identity-preserving exports and audited DP4+/STS formula implementations."""
import csv
import json
import math
from pathlib import Path
import tempfile
import time

from nmrflow import ROOT, require, save, digest, gaussian_read, rms, calculation_config


def table(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def write_csv(path, columns, rows):
    with Path(path).open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=columns)
        w.writeheader()
        w.writerows(rows)


def beta_fraction(a, b, x):
    """Lentz continued fraction for the regularized incomplete beta."""
    qab, qap, qam = a+b, a+1, a-1
    c, d = 1.0, 1-qab*x/qap
    if abs(d) < 1e-300:
        d = 1e-300
    d, h = 1/d, 1/d
    for m in range(1, 10001):
        aa = m*(b-m)*x/((qam+2*m)*(a+2*m))
        d, c = 1+aa*d, 1+aa/c
        if abs(d) < 1e-300:
            d = 1e-300
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1/d
        h *= d*c
        aa = -(a+m)*(qab+m)*x/((a+2*m)*(qap+2*m))
        d, c = 1+aa*d, 1+aa/c
        if abs(d) < 1e-300:
            d = 1e-300
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1/d
        factor = d*c
        h *= factor
        if abs(factor-1) < 3e-14:
            return h
    raise ValueError('Incomplete beta did not converge')


def student_tail(x, df):
    """One-sided Student-t survival at abs(x), equivalent to Excel TDIST(x,df,1)."""
    require(df > 0 and math.isfinite(x), 'Invalid t distribution argument')
    x = abs(x)
    if x == 0:
        return .5
    z = df/(df+x*x)
    a, b = df/2, .5
    if z == 0:
        return 0.0
    if z == 1:
        return .5
    factor = math.exp(math.lgamma(a+b)-math.lgamma(a)-math.lgamma(b)+a*math.log(z)+b*math.log1p(-z))
    if z < (a+1)/(a+b+2):
        val = factor*beta_fraction(a,b,z)/a
    else:
        val = 1-factor*beta_fraction(b,a,1-z)/b
    return max(0.0, min(.5, val/2))


def linear_fit(x, y):
    require(len(x) == len(y) and len(x) >= 3, 'At least three assigned signals per fitted nucleus required')
    mx, my = sum(x)/len(x), sum(y)/len(y)
    denominator = sum((v-mx)**2 for v in x)
    require(denominator > 1e-15, 'Degenerate regression')
    slope = sum((a-mx)*(b-my) for a,b in zip(x,y))/denominator
    require(abs(slope) > 1e-12, 'Zero regression slope')
    return slope, my-slope*mx


def normalize(log_values):
    m = max(log_values.values())
    require(math.isfinite(m), 'All likelihoods underflow; parameter/data check required')
    factors = {k: math.exp(v-m) for k,v in log_values.items()}
    total = sum(factors.values())
    return {k:100*v/total for k,v in factors.items()}


def dp4plus(data, params):
    """Workbook convention: regress calculated shifts on experimental, then invert."""
    log_scores, detail = {}, {}
    for cid, rows in data.items():
        scores = {'H_scaled':0.0,'H_unscaled':0.0,'C_scaled':0.0,'C_unscaled':0.0}
        used = set()
        fit_info = {}
        for element in ('H','C'):
            subset = [r for r in rows if r['nucleus'] == element]
            if not subset:
                continue
            used.add(element)
            exp = [r['experimental_ppm'] for r in subset]
            unscaled = [params['TMS'][element]-r['shielding'] for r in subset]
            slope, intercept = linear_fit(exp, unscaled)
            fit_info[element] = {'slope':slope,'intercept':intercept}
            for row, u in zip(subset, unscaled):
                hybrid = row['hybridization']
                require(hybrid in ('sp2','sp3'), 'DP4+ requires sp2/sp3 classification for every C/H signal')
                mu, sigma, df = params['unscaled'][element+'_'+hybrid]
                s_sigma, s_df = params['scaled'][element]
                error_u = u-row['experimental_ppm']
                error_s = (u-intercept)/slope-row['experimental_ppm']
                # Excel's legacy TDIST truncates degrees of freedom to integers.
                pu = student_tail((error_u-mu)/sigma,int(df))
                ps = student_tail(error_s/s_sigma,int(s_df))
                require(pu > 0 and ps > 0, 'Probability underflow; inspect extreme errors')
                scores[element+'_unscaled'] += math.log(pu)
                scores[element+'_scaled'] += math.log(ps)
        require(bool(used), 'No DP4+ signals')
        scores = {k:v for k,v in scores.items() if k[0] in used}
        for element in used:
            scores[element+'_combined'] = scores[element+'_scaled']+scores[element+'_unscaled']
        scores['all_combined'] = sum(scores[e+'_combined'] for e in used)
        log_scores[cid] = scores
        detail[cid] = {'fits':fit_info,'log_likelihoods':scores}
    keys = set.intersection(*(set(v) for v in log_scores.values()))
    return {'probabilities_percent':{k:normalize({c:v[k] for c,v in log_scores.items()}) for k in sorted(keys)},
            'details':detail}


def dp4(data, params):
    """Classic Goodman DP4: internally scaled errors and Gaussian two-tailed probabilities.

    Fitting raw shieldings and inverting is algebraically equivalent to fitting
    TMS-referenced shifts and inverting; the reference cancels from this model.
    """
    scores, fits = {}, {}
    for cid,rows in data.items():
        scores[cid],fits[cid] = {},{}
        for nucleus in ['C','H']:
            subset = [r for r in rows if r['nucleus']==nucleus]
            if not subset:
                continue
            slope,intercept = linear_fit([r['experimental_ppm'] for r in subset],
                                          [r['shielding'] for r in subset])
            logs = []
            for row in subset:
                error = (row['shielding']-intercept)/slope-row['experimental_ppm']
                z = abs((error-params['mean'][nucleus])/params['sigma'][nucleus])
                probability = math.erfc(z/math.sqrt(2))
                require(probability > 0,'Classic DP4 normal tail underflow; inspect errors')
                logs.append(math.log(probability))
            scores[cid][nucleus] = sum(logs)
            fits[cid][nucleus] = {'slope':slope,'intercept':intercept}
        require(bool(scores[cid]),'No classic DP4 signals')
        scores[cid]['combined'] = sum(scores[cid].values())
    keys = set.intersection(*(set(x) for x in scores.values()))
    return {'probabilities_percent':{k:normalize({c:v[k] for c,v in scores.items()}) for k in sorted(keys)},
            'fits':fits,'model':'Classic DP4, Goodman normal-error model'}


def sts(data, params):
    result, logp = {}, {}
    for cid, rows in data.items():
        subset = [r for r in rows if r['nucleus'] == 'C']
        require(len(subset) >= 3, 'STS requires at least three carbon signals')
        predictions = []
        for row in subset:
            kind, heavy = row['carbon_type'], row['heavy_atom'] or 'None'
            require(kind in params['scale'], 'Unknown STS carbon type '+kind)
            value = row['shielding']
            if heavy != 'None':
                key = kind+'|'+heavy
                require(key in params['heavy'], 'Unsupported STS heavy-atom correction '+key)
                a,b = params['heavy'][key]
                value = a*value+b
            a,b = params['scale'][kind]
            predictions.append(a*value+b)
        exp = [r['experimental_ppm'] for r in subset]
        # STS workbook fits experimental shifts on the segment-corrected predictions.
        slope, intercept = linear_fit(predictions,exp)
        errors = [slope*p+intercept-e for p,e in zip(predictions,exp)]
        logs = []
        for row,error in zip(subset,errors):
            sigma,df = params['error'][row['carbon_type']]
            p = 2*student_tail(error/sigma,df)
            require(p > 0, 'STS tail underflow')
            logs.append(math.log(p))
        logp[cid] = sum(logs)
        result[cid] = {'MAE':sum(abs(e) for e in errors)/len(errors),
                       'RMSE':math.sqrt(sum(e*e for e in errors)/len(errors)),
                       'Pmean':math.exp(sum(logs)/len(logs)), 'slope':slope,'intercept':intercept}
    relative = normalize(logp)
    for cid in result:
        result[cid]['Prel_percent'] = relative[cid]
    return result


def export(cfg, output, config_dir):
    output.mkdir(parents=True, exist_ok=True)
    exports = Path(tempfile.mkdtemp(prefix='.exports_pending_',dir=str(output)))
    averaged, raw, provenance = {}, [], []
    for candidate in cfg['candidates']:
        cid = candidate['id']
        manifest = output/cid/'conformers.json'
        state = json.loads(manifest.read_text())
        require(calculation_config(state['config']) == calculation_config(cfg), 'Manifest/config mismatch for '+cid)
        ids = [r['id'] for r in state['records']]
        require(len(ids) == len(set(ids)) and len(state['nmr_ids']) == len(set(state['nmr_ids']))
                and set(state['nmr_ids']) <= set(ids), 'Missing/duplicate conformer ID')
        chosen = [r for r in state['records'] if r['id'] in state['nmr_ids']]
        require(all(math.isfinite(r['nmr_weight']) and r['nmr_weight'] >= 0 for r in chosen), 'Invalid NMR weights')
        require(chosen and abs(sum(r['nmr_weight'] for r in chosen)-1) < 1e-8, 'Invalid NMR weight sum')
        averaged[cid] = {}
        for method in cfg['methods']:
            sums, identity = {}, None
            for row in chosen:
                parsed = gaussian_read(row['nmr'][method], 'nmr')
                receipt = Path(row['nmr'][method]).parent/'receipt.json'
                require(receipt.exists() and json.loads(receipt.read_text())['output_sha256'] == digest(row['nmr'][method]),
                        'NMR stage receipt/hash mismatch')
                require(rms(parsed['geometry'],row['geometry']) < 1e-4, 'NMR geometry/manifest mismatch')
                current = {i:el for i,(el,_) in parsed['shielding'].items()}
                if identity is None:
                    identity = current
                require(current == identity, 'NMR atom set mismatch')
                source_hash = digest(row['nmr'][method])
                for atom,(element,shield) in parsed['shielding'].items():
                    sums[atom] = sums.get(atom,0.0)+row['nmr_weight']*shield
                    raw.append({'candidate':cid,'method':method,'conformer':row['id'],'atom':atom,'element':element,
                                'shielding_ppm':shield,'weight':row['nmr_weight'],'source_sha256':source_hash})
            averaged[cid][method] = {a:(identity[a],s) for a,s in sums.items()}
            # Same atom-number / element / weighted raw shielding layout used by upstream result files.
            contributions = '\n'.join('contribution of conformer {}: {:.10f}'.format(i,r['nmr_weight'])
                                      for i,r in enumerate(chosen,1))+'\n'
            (exports/(cid+'_nmr'+method+'_result.txt')).write_text(contributions+'\n'.join(
                '{}\t{}\t{:.10f}'.format(a,identity[a],sums[a]) for a in sorted(sums) if identity[a] in ('C','H'))+'\n')
            (exports/(cid+'_nmr'+method+'_carbon_result.txt')).write_text(contributions+'\n'.join(
                '{}\tC\t{:.10f}'.format(a,sums[a]) for a in sorted(sums) if identity[a]=='C')+'\n')
        provenance.append({'candidate':cid,'manifest_sha256':digest(manifest),
                           'coverage_within_computed_ensemble':state['coverage_within_computed_ensemble'],
                           'ensemble_total':state['ensemble_total'],'dft_selected':state['dft_selected']})
    write_csv(exports/'per_conformer_shielding.csv',list(raw[0]),raw)
    average_rows = [{'candidate':c,'method':m,'atom':a,'element':el,'shielding_ppm':s}
                    for c,methods in averaged.items() for m,atoms in methods.items() for a,(el,s) in atoms.items()]
    write_csv(exports/'weighted_shielding.csv',list(average_rows[0]),average_rows)
    assignments = config_dir/cfg['assignments']
    experiments = config_dir/cfg['experiment']
    if not assignments.exists() or not experiments.exists():
        save(exports/'comparison_status.json', {'status':'raw_exports_complete_assignment_required',
             'needed':[str(assignments),str(experiments)],'provenance':provenance})
        print('Raw exports complete. Experimental assignments required for DP4+/STS.')
        publish(exports,output)
        return
    exp = table(experiments)
    require(exp and len(set(r['signal'] for r in exp)) == len(exp), 'Empty/duplicate experimental signal ID')
    require(all(r['nucleus'] in ('C','H') for r in exp), 'Only C/H experimental signals supported')
    mapping = table(assignments)
    require(len(set((r['candidate'],r['signal']) for r in mapping)) == len(mapping), 'Duplicate atom assignment')
    lookup = {(r['candidate'],r['signal']):r for r in mapping}
    assembled = {}
    for method in cfg['methods']:
        data = {}
        for candidate in cfg['candidates']:
            cid = candidate['id']
            data[cid], used = [], set()
            for e in exp:
                require((cid,e['signal']) in lookup, 'Missing experimental assignment '+cid+':'+e['signal'])
                row = lookup[(cid,e['signal'])]
                atoms = [int(x) for x in row['atoms'].split(';')]
                require(atoms and len(set(atoms)) == len(atoms), 'Invalid equivalence group')
                require(not (used & set(atoms)), 'Atom belongs to multiple experimental signals')
                used.update(atoms)
                require(all(a in averaged[cid][method] and averaged[cid][method][a][0] == e['nucleus'] for a in atoms), 'Assignment element/index mismatch')
                ppm = float(e['ppm'])
                require(math.isfinite(ppm), 'Nonfinite experimental shift')
                data[cid].append(dict(row, nucleus=e['nucleus'],experimental_ppm=ppm,
                                     shielding=sum(averaged[cid][method][a][1] for a in atoms)/len(atoms)))
        assembled[method] = data
        rows = [{'signal':e['signal'],'nucleus':e['nucleus'],'experimental_ppm':float(e['ppm']),
                 **{c:next(r['shielding'] for r in d if r['signal']==e['signal']) for c,d in data.items()}} for e in exp]
        write_csv(exports/('method'+method+'_comparison_input.csv'),list(rows[0]),rows)
    parameters = json.loads((ROOT/'calibration.json').read_text())
    descriptive = {}
    for method,data in assembled.items():
        descriptive[method] = {}
        for cid,rows in data.items():
            descriptive[method][cid] = {}
            for nucleus in ['C','H']:
                subset = [r for r in rows if r['nucleus']==nucleus]
                if len(subset) < 3:
                    continue
                a,b = linear_fit([r['shielding'] for r in subset],[r['experimental_ppm'] for r in subset])
                errors = [a*r['shielding']+b-r['experimental_ppm'] for r in subset]
                descriptive[method][cid][nucleus] = {'slope':a,'intercept':b,
                    'fitted_MAE':sum(abs(e) for e in errors)/len(errors),
                    'fitted_RMSE':math.sqrt(sum(e*e for e in errors)/len(errors)),
                    'note':'Fit to this experimental dataset; descriptive only, not an external validation.'}
    statistics = {}
    for branch in cfg['postprocessing']:
        method, kind = branch['method'], branch['kind']
        require(method in assembled, 'Postprocessing branch needs disabled method '+method)
        validate_branch_data(branch,assembled[method],parameters)
        # Selected geometry/weighting and optional solvent/defaults differ from calibration.
        compatibility = 'exploratory_custom_protocol'
        if not branch.get('allow_exploratory',False):
            statistics[kind+'_'+method] = {'status':'input_exported_calibration_review_required',
                'parameter_key':branch['parameter_key'],'differences':branch['differences']}
            continue
        if kind in ('dp4','dp4plus'):
            require(len(assembled[method]) >= 2, 'DP4+ relative comparison requires at least two candidates')
            if kind == 'dp4':
                params = parameters['dp4_classic'][branch['parameter_key']]
                result = dp4(assembled[method],params)
            else:
                params = parameters['dp4'][branch['parameter_key']]
                result = dp4plus(assembled[method],params)
        elif kind == 'sts':
            params = parameters['sts'][branch['parameter_key']]
            result = sts(assembled[method],params)
        else:
            raise ValueError('Unknown postprocessor '+kind)
        statistics[kind+'_'+method] = {'status':compatibility,'parameter_key':branch['parameter_key'],
             'differences':branch['differences'],'result':result}
    save(exports/'comparison_results.json', {'protocol':cfg['protocol'],'statistics':statistics,
         'method_layout':cfg.get('method_layout','historical_v2_v3'),
         'nmr_settings':cfg.get('nmr_settings',{'solvent_model':'SMD','gaussian_defaults':'G16'}),
         'nmr_profile':cfg.get('nmr_profile','historical_smd_g16'),
         'clustering_profile':cfg.get('clustering_profile','historical_fixed'),
         'clustering':cfg['clustering'],
         'descriptive_fits':descriptive,
         'provenance':provenance,'calibration_sha256':digest(ROOT/'calibration.json'),
         'assignments_sha256':digest(assignments),'experiment_sha256':digest(experiments),
         'interpretation':'Relative comparisons within the supplied candidates; custom protocol is not a validated standard probability.'})
    workbook_exports(cfg,assembled,exports)
    publish(exports,output)
    print('NMR exports and comparisons complete: '+str(output/'exports'))


def publish(staging, output):
    final = output/'exports'
    if final.exists():
        final.rename(output/('exports_previous_'+time.strftime('%Y%m%d_%H%M%S')+'_'+staging.name[-6:]))
    staging.rename(final)


def validate_branch_data(branch, data, parameters):
    """Validate classifications even when native exploratory scores are disabled."""
    kind = branch['kind']
    for rows in data.values():
        if kind == 'dp4plus':
            require(all(r['hybridization'] in ('sp2','sp3') for r in rows), 'DP4+ requires explicit sp2/sp3 classification')
            for nucleus in ('C','H'):
                subset = [r for r in rows if r['nucleus'] == nucleus]
                require(not subset or len(subset) >= 3, 'DP4+ requires at least three assigned signals per fitted nucleus')
        elif kind == 'sts':
            params = parameters['sts'][branch['parameter_key']]
            carbon = [r for r in rows if r['nucleus'] == 'C']
            require(len(carbon) >= 3, 'STS requires at least three carbon signals')
            for r in carbon:
                k,h = r['carbon_type'],r['heavy_atom'] or 'None'
                require(k in params['scale'], 'Unknown STS carbon type '+k)
                require(h == 'None' or k+'|'+h in params['heavy'], 'Unsupported STS heavy-atom correction '+k+'|'+h)


def workbook_exports(cfg, assembled, exports):
    import openpyxl
    for branch in cfg['postprocessing']:
        data = assembled[branch['method']]
        candidates = list(data)
        first = data[candidates[0]]
        if branch['kind'] == 'sts':
            book = openpyxl.load_workbook(ROOT/'upstream/STS.xlsx')
            sheet = book[branch['parameter_key']]
            carbon = [r for r in first if r['nucleus'] == 'C']
            require(len(carbon) <= 91 and len(candidates) <= 12, 'STS workbook capacity exceeded')
            for n,row in enumerate(carbon,9):
                sheet.cell(n,3,row['heavy_atom'] or 'None')
                sheet.cell(n,4,row['carbon_type'])
                sheet.cell(n,5,row['signal'])
                sheet.cell(n,6,row['experimental_ppm'])
                for j,c in enumerate(candidates):
                    other = next(x for x in data[c] if x['signal']==row['signal'])
                    require(other['carbon_type']==row['carbon_type'] and other['heavy_atom']==row['heavy_atom'], 'STS workbook shares carbon classes across candidates; classify consistently')
                    sheet.cell(n,7+6*j,other['shielding'])
                    sheet.cell(7,7+6*j,c)
        elif branch['kind'] == 'dp4plus':
            book = openpyxl.load_workbook(ROOT/'upstream/DP4plus.xlsx')
            sheet = book['Main']
            params = json.loads((ROOT/'calibration.json').read_text())['dp4'][branch['parameter_key']]
            sheet['A6'],sheet['C6'],sheet['E6'] = params['functional'],params['solvation'],params['basis']
            sheet['A2'],sheet['D2'] = 'Default','Shielding tensors'
            # The source workbook contains one hardcoded J6; restore table-based value.
            sheet['J6'] = params['unscaled']['H_sp3'][1]
            require(len(first)<=300 and len(candidates)<=128,'DP4+ workbook capacity exceeded')
            for n,row in enumerate(first,14):
                sheet.cell(n,1,row['nucleus'])
                sheet.cell(n,2,'x' if row['hybridization']=='sp2' else None)
                sheet.cell(n,3,row['experimental_ppm'])
                for j,c in enumerate(candidates):
                    other = next(x for x in data[c] if x['signal']==row['signal'])
                    require(other['hybridization']==row['hybridization'],'DP4 workbook shares hybridization across candidates')
                    sheet.cell(n,4+j,other['shielding'])
            for j,c in enumerate(candidates):
                sheet.cell(9,4+j,c)
                sheet.cell(13,4+j,c)
                book['Detailed Results'].cell(2,3+j,c)
        else:
            require(branch['kind']=='dp4','Unknown workbook branch')
            book = openpyxl.Workbook()
            sheet = book.active
            sheet.title = 'DP4_input'
            sheet.append(['Signal','Nucleus','Experimental ppm']+candidates)
            for row in first:
                sheet.append([row['signal'],row['nucleus'],row['experimental_ppm']]+
                             [next(r['shielding'] for r in data[c] if r['signal']==row['signal']) for c in candidates])
            sheet.cell(1,len(candidates)+4,'Candidate columns contain weighted raw shieldings (ppm)')
        notes = book.create_sheet('Protocol_metadata')
        notes.append(['Status','exploratory_custom_protocol'])
        notes.append(['Method',branch['method']])
        from protocol_options import method_map
        functional,basis,kind=method_map(cfg)[branch['method']]
        notes.append(['NMR level',functional+'/'+basis])
        notes.append(['Method layout',cfg.get('method_layout','historical_v2_v3')])
        notes.append(['NMR settings',json.dumps(cfg.get('nmr_settings',{'solvent_model':'SMD','gaussian_defaults':'G16'}))])
        notes.append(['Clustering',json.dumps(cfg['clustering'])])
        notes.append(['Differences','; '.join(branch['differences'])])
        notes.append(['Note','Open in Excel to recalculate workbook formulas; JSON results are computed independently.'])
        book.calculation.fullCalcOnLoad = True
        book.save(exports/(branch['kind']+'_method'+branch['method']+'_exploratory.xlsx'))
