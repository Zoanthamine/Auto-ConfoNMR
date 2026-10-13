"""Cross-version checks: V1.2 batch interfaces with V2 optional modules."""
import configparser
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import experiment_batch as batch
import module_options as mo
import nmrflow as nf
import project_console as pc
import progress_status as ps


class IntegrationChecks(unittest.TestCase):
    def test_module_csvs_and_templates_are_not_nmr_samples(self):
        cfg=json.loads((nf.ROOT/'job.example.json').read_text())
        cfg['module_analysis']=dict(mo.DEFAULT_ANALYSIS,j_experiment='experiment_customj.csv')
        settings=configparser.ConfigParser();settings.read_string(pc.settings_text(cfg))
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            for name in ['experiment_j.csv','experiment_ecd.csv','experiment_customj.csv','experiment_j.template.csv','experiment_ecd.template.csv']:
                (root/name).write_text('different,module,schema\n')
            self.assertEqual(batch.discover(root,settings,allow_empty=True),[])
            (root/'assignments.csv').write_text('mapping\n')
            (root/'experiment_product.csv').write_text('signal,nucleus,ppm\n')
            self.assertEqual([r['sample'] for r in batch.discover(root,settings)],['product'])

    def test_standalone_check_only_does_not_export_or_submit(self):
        cfg=json.loads((nf.ROOT/'job.example.json').read_text())
        cfg['modules']=dict(mo.DEFAULT_MODULES,base='conformers',coupling=True)
        p=configparser.ConfigParser();p.read_string(pc.settings_text(cfg))
        with tempfile.TemporaryDirectory() as d, patch.object(pc,'sync',return_value=(Path(d)/'job.json',cfg,p)), patch.object(pc.subprocess,'run') as run, patch.object(pc,'refresh_modules') as refresh, patch.object(pc,'queue_guard') as queue, contextlib.redirect_stdout(io.StringIO()):
            result=pc.compare(Path(d),check_only=True)
            self.assertEqual(result['status'],'no_nmr_experiment')
            run.assert_not_called();refresh.assert_not_called();queue.assert_not_called()

    def test_live_coupling_counts_without_fake_nmr(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);cfg={'output':'raw','methods':['01','02','03'],'candidates':[{'id':'001RR'}],
                 'modules':dict(mo.DEFAULT_MODULES,base='conformers',coupling=True)}
            for stage in ['coupling','coupling_shielding']:
                folder=root/'raw/001RR/conformers/c0001'/stage;folder.mkdir(parents=True)
                (folder/'receipt.json').write_text(json.dumps({'status':'accepted'}))
            snap=ps.snapshot(root,cfg);row=snap['candidates'][0]
            self.assertEqual(row['nmr'],{})
            self.assertEqual(row['properties'],{'coupling':1,'coupling_shielding':1})
            self.assertEqual(row['property_complete_conformers'],1)
            text=ps.format_snapshot(snap)
            self.assertNotIn('NMR任务',text);self.assertIn('coupling=1',text)

    def test_quick_progress_identifies_weight_source(self):
        with tempfile.TemporaryDirectory() as d:
            cfg={'output':'raw','candidates':[{'id':'001RR'}],
                 'modules':dict(mo.DEFAULT_MODULES,base='conformers',conformer_weight='gfn2_energy')}
            snap=ps.snapshot(Path(d),cfg)
            self.assertIn('GFN2电子能估算',ps.format_snapshot(snap))
            self.assertEqual(snap['candidates'][0]['properties'],{})


if __name__=='__main__':unittest.main()
