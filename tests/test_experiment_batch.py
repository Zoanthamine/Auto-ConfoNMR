"""Multi-sample input checks, isolated comparison exports and complete collection."""
import configparser
import json
from pathlib import Path
import re
import sys
import tempfile
import unittest
from unittest.mock import patch

import experiment_batch as batch
import nmrflow as nf
import postprocess as pp
import project_console as pc
from test_nmrflow import gaussian


class SampleChecks(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.cfg = json.loads((nf.ROOT/'job.example.json').read_text(encoding='utf-8'))
        self.cfg.update(candidates=[{'id': 'candidate_A', 'input': 'A.xyz'}, {'id': 'candidate_B', 'input': 'B.xyz'}], output='raw')
        self.settings = configparser.ConfigParser(interpolation=None)
        self.settings.read_string(pc.settings_text(self.cfg))
        self.out = self.root/'raw'
        self.geom = [['C',0,0,0],['C',1.5,0,0],['C',3,0,0],['H',0,1,0],['H',0,0,1],['H',0,-1,0],['H',1.5,1,0],['H',3,1,0]]
        for candidate in self.cfg['candidates']:
            folder = self.out/candidate['id']/'sampling'; folder.mkdir(parents=True)
            nf.save(folder/'atom_map.json', [{'calculation_atom': i, 'input_serial': 10+i, 'element': atom[0]} for i, atom in enumerate(self.geom, 1)])
        (self.root/'assignments.csv').write_bytes((nf.ROOT/'assignments.example.csv').read_bytes())
        for name in ['A', 'B', 'C']:
            (self.root/('experiment_'+name+'.csv')).write_bytes((nf.ROOT/'experiment.example.csv').read_bytes())
        self.items = batch.discover(self.root, self.settings)

    def tearDown(self):
        self.tmp.cleanup()

    def test_unlimited_discovery_specific_mapping_and_single(self):
        self.assertEqual([i['sample'] for i in self.items], ['A', 'B', 'C'])
        (self.root/'assignments_B.csv').write_bytes((self.root/'assignments.csv').read_bytes())
        (self.root/'experiment.csv').write_bytes((self.root/'experiment_A.csv').read_bytes())
        items = batch.discover(self.root, self.settings)
        self.assertEqual(len(items), 4)
        self.assertEqual(next(i for i in items if i['sample']=='B')['assignments'].name, 'assignments_B.csv')
        self.assertEqual(next(i for i in items if i['sample']=='C')['assignments'].name, 'assignments.csv')
        for n in range(20):
            (self.root/('experiment_product_'+str(n)+'.csv')).write_bytes((self.root/'experiment_A.csv').read_bytes())
        self.assertEqual(len(batch.discover(self.root, self.settings)), 24)

    def test_legacy_override_and_name_collision(self):
        self.settings.set('对比', 'upper_experiment', 'experiment_A.csv')
        self.settings.set('对比', 'upper_assignments', 'assignments.csv')
        items = batch.discover(self.root, self.settings)
        self.assertEqual(len(items), 3)
        self.assertIn('upper', [i['sample'] for i in items])
        (self.root/'experiment_UPPER.csv').write_bytes((self.root/'experiment_A.csv').read_bytes())
        with self.assertRaisesRegex(ValueError, '冲突'):
            batch.discover(self.root, self.settings)

    def test_bad_header_nonfinite_and_duplicate_signal(self):
        exp = self.items[0]['experiment']; original = exp.read_text(encoding='utf-8-sig')
        for text, error in [(original.replace('ppm','PPM'), '表头'),
                            (original.replace('15.2',' 15.2'), '首尾空格'),
                            (original.replace('15.2','NaN'), 'NaN/Inf'),
                            (original+'C1,C,12\n', '重复')]:
            exp.write_text(text, encoding='utf-8')
            with self.assertRaisesRegex(ValueError, error):
                batch.validate_sample(self.items[0], self.cfg, self.out)
        exp.write_text(original, encoding='utf-8')

    def test_missing_extra_overlap_wrong_nucleus_and_classes(self):
        path = self.root/'assignments.csv'; original = path.read_text(encoding='utf-8-sig')
        for text, error in [(original.replace('candidate_B,H3,8,sp2,,\n',''), '缺少'),
                            (original+'stray,C1,1,sp3,sp3-CH3,None\n', '多余'),
                            (original.replace('H1,4;5;6','H1,4;4;6'), '重复'),
                            (original.replace('H1,4;5;6','H1,9'), '核种'),
                            (original.replace('candidate_B,C1,1,sp3','candidate_B,C1,1,sp2'), '分类不一致'),
                            (original.replace('sp3-CH3','bad-type'), 'Unknown STS')]:
            path.write_text(text, encoding='utf-8')
            with self.assertRaisesRegex(ValueError, error):
                batch.validate_sample(self.items[0], self.cfg, self.out)
        path.write_text(original, encoding='utf-8')

    def test_equivalence_manual_review_and_fail_all_preflight(self):
        report = batch.validate_sample(self.items[0], self.cfg, self.out)
        self.assertTrue(any('等价原子' in x for x in report['warnings']))
        (self.root/'experiment_B.csv').write_text('signal,nucleus,ppm\nC1,C,NaN\n', encoding='utf-8')
        with patch.object(pp, 'export') as export:
            with self.assertRaisesRegex(ValueError, '未进行本批比较'):
                batch.run_batch(self.root, self.cfg, self.settings)
            export.assert_not_called()
        self.assertFalse((self.root/'results').exists())

    def test_check_only_and_incomplete_do_not_export(self):
        with patch.object(pp, 'export') as export:
            result = batch.run_batch(self.root, self.cfg, self.settings, check_only=True)
            self.assertEqual(len(result['samples']), 3)
            with self.assertRaisesRegex(ValueError, '未完成'):
                batch.run_batch(self.root, self.cfg, self.settings)
            export.assert_not_called()
        self.assertFalse((self.out/'.workflow.guard').exists())

    def test_entries_forward_check_flag_and_ini(self):
        pc.install_entries(self.root, nf.ROOT, '/migration/env/bin/python')
        shell = (self.root/'03_compare.sh').read_text(encoding='utf-8')
        self.assertIn('"$@"', shell)
        self.assertIn('--check-only', shell)
        self.assertIn('experiment_pattern = experiment_*.csv', pc.settings_text(self.cfg))

    def completed_fixture(self):
        cfg = json.loads(json.dumps(self.cfg)); cfg['output'] = str(self.out)
        for c in cfg['candidates']:
            c['input'] = str(self.root/c['input'])
        for c in cfg['candidates']:
            logs = {}
            for method in cfg['methods']:
                work = self.out/c['id']/'conformers/c0001'/('nmr'+method); work.mkdir(parents=True)
                text = gaussian(self.geom)
                if c['id']=='candidate_B':
                    text = re.sub(r'Isotropic = ([-\d.]+)', lambda m: 'Isotropic = '+str(float(m[1])+0.3), text)
                log = work/'job.log'; log.write_text(text)
                nf.save(work/'receipt.json', {'output_sha256': nf.digest(log)})
                logs[method] = str(log)
            nf.save(self.out/c['id']/'conformers.json', {'config': cfg, 'records': [dict(id='c0001', geometry=self.geom, weight=1, gibbs_hartree=-40.0, nmr_weight=1, nmr=logs)],
                    'nmr_ids': ['c0001'], 'coverage_within_computed_ensemble': 1, 'ensemble_total': 1, 'dft_selected': 1})
            nf.save(self.out/c['id']/'state.json', {'status': 'complete', 'identity': {'config': nf.calculation_config(cfg)}})
        return cfg

    @unittest.skipUnless(all((nf.ROOT/'upstream'/name).is_file() for name in ('DP4plus.xlsx', 'STS.xlsx')), 'Upstream resource workbooks required')
    def test_three_samples_nine_workbooks_isolated_and_collected(self):
        self.completed_fixture()
        (self.root/'experiment_B.csv').write_text((self.root/'experiment_B.csv').read_text(encoding='utf-8-sig').replace('15.2','17.2'), encoding='utf-8')
        old = self.out/'exports'; old.mkdir(); (old/'old.txt').write_text('retained')
        (self.root/'job.json').write_text(json.dumps(self.cfg),encoding='utf-8')
        (self.root/'project.ini').write_text(pc.settings_text(self.cfg),encoding='utf-8')
        with patch.object(pc,'queue_guard'):
            result = pc.compare(self.root)
        self.assertEqual(result['status'], 'complete')
        self.assertEqual(len(result['samples']), 3)
        self.assertEqual((old/'old.txt').read_text(encoding='utf-8'), 'retained')
        hashes = []
        for entry in result['samples']:
            folder = self.root/entry['result']
            self.assertEqual(len(list(folder.glob('*.xlsx'))), 3)
            data = json.loads((folder/'comparison_results.json').read_text(encoding='utf-8'))
            self.assertEqual(set(data['statistics']), {'dp4plus_01','dp4plus_02','sts_03'})
            hashes.append(data['experiment_sha256'])
            self.assertTrue((folder/'数据与结果复核.md').is_file())
            self.assertTrue((folder/'信号误差.csv').is_file())
        self.assertNotEqual(hashes[0], hashes[1]); self.assertEqual(hashes[0], hashes[2])
        rows = batch.collect_comparisons(self.root)
        self.assertEqual(len(rows), 3)
        self.assertEqual(len(list((self.root/'最终结果/实验比较').rglob('*.xlsx'))), 9)
        self.assertTrue(all(not f.is_symlink() for f in (self.root/'最终结果/实验比较').rglob('*')))
        self.assertEqual(len(batch.collect_comparisons(self.root)), 3)
        self.assertEqual(len(pp.table(self.root/'最终结果/各样品最新结果.csv')), 3)
        self.assertEqual(len(list((self.root/'最终结果/实验比较').rglob('输入记录/experiment.csv'))), 3)

    def test_partial_failure_keeps_success_and_history(self):
        self.completed_fixture()
        calls = []
        def fake_export(cfg, out, config_dir, destination=None):
            calls.append(destination.name)
            if destination.name == 'B':
                raise ValueError('diagnostic failure')
            destination.mkdir()
            (destination/'dummy.txt').write_text('A completed')
        with patch.object(pp, 'export', side_effect=fake_export), patch.object(batch, 'review_results'):
            with self.assertRaisesRegex(ValueError, 'diagnostic failure'):
                batch.run_batch(self.root, self.cfg, self.settings)
        self.assertEqual(calls, ['A','B'])
        manifest = next((self.root/'results').glob('comparison_*/batch.json'))
        status = json.loads(manifest.read_text(encoding='utf-8'))
        self.assertEqual(status['status'], 'failed')
        self.assertEqual([e['status'] for e in status['samples']], ['complete','failed'])
        self.assertFalse((self.root/'results/latest_batch.json').exists())
        self.assertEqual(len(batch.collect_comparisons(self.root)), 1)
        # A later run has a distinct directory and retains the failed run's evidence.
        with patch.object(pp, 'export', side_effect=fake_export), patch.object(batch, 'review_results'):
            with self.assertRaises(ValueError):
                batch.run_batch(self.root, self.cfg, self.settings)
        self.assertEqual(len(list((self.root/'results').glob('comparison_*/batch.json'))), 2)
        self.assertEqual(len(batch.collect_comparisons(self.root)), 2)


if __name__ == '__main__':
    unittest.main()
