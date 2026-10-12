"""Regression checks for log preservation, configuration identity and queue safety."""
import contextlib
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import nmrflow as nf
import project_console as pc
import project_views as pv


class ReadabilityChecks(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.cfg=json.loads((nf.ROOT/'job.example.json').read_text())
        self.cfg['candidates']=[{'id':'A','input':'inputs/A.pdb'}];self.cfg['output']='raw'
        (self.root/'inputs').mkdir();(self.root/'inputs/A.pdb').write_text('COMPND A\n')
        (self.root/'job.json').write_text(json.dumps(self.cfg))
        (self.root/'项目配置.ini').write_text(pc.settings_text(self.cfg),encoding='utf-8')

    def tearDown(self):self.tmp.cleanup()

    def test_stream_preserves_software_bytes_and_batch_stdin(self):
        log=self.root/'software.log';capture=io.StringIO()
        with patch.dict(os.environ,{'NMR_STREAM_LOG':'1','NMR_PROJECT_DIR':''}),contextlib.redirect_stdout(capture):
            nf.command([sys.executable,'-u','-c',"import sys;print('INPUT='+sys.stdin.read().strip());print('ERR',file=sys.stderr)"],self.root,log,stdin='menu\n')
        self.assertEqual(log.read_text(),'INPUT=menu\nERR\n')
        self.assertIn('INPUT=menu',capture.getvalue());self.assertIn('软件结束',capture.getvalue())

    def test_failure_keeps_output_and_returncode(self):
        log=self.root/'failed.log'
        with patch.dict(os.environ,{'NMR_STREAM_LOG':'1','NMR_PROJECT_DIR':''}),contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaisesRegex(ValueError,'Command failed'):
                nf.command([sys.executable,'-u','-c',"import sys;print('DIAGNOSIS');sys.exit(7)"],self.root,log)
        self.assertEqual(log.read_text(),'DIAGNOSIS\n')

    def test_chinese_settings_roundtrip(self):
        _,cfg,_=pc.sync(self.root)
        self.assertEqual(cfg['max_conformers'],6);self.assertTrue(all(x['allow_exploratory'] for x in cfg['postprocessing']))
        self.assertEqual(cfg['shermo']['T'],298.15)

    def test_existing_computation_rejects_parameter_change_without_overwrite(self):
        _,cfg,_=pc.sync(self.root);resolved=json.loads(json.dumps(cfg));resolved['output']=str(self.root/'raw')
        resolved['candidates'][0]['input']=str(self.root/'inputs/A.pdb')
        state=self.root/'raw/A/state.json';state.parent.mkdir(parents=True)
        state.write_text(json.dumps({'identity':{'config':nf.calculation_config(resolved)},'status':'complete'}))
        pc.sync(self.root)  # An unchanged calculation may change only analysis settings.
        previous=(self.root/'job.json').read_bytes()
        ini=self.root/'项目配置.ini';ini.write_text(ini.read_text(encoding='utf-8').replace('max_conformers = 6','max_conformers = 7'),encoding='utf-8')
        with self.assertRaisesRegex(ValueError,'计算配置发生变化'):pc.sync(self.root)
        self.assertEqual(previous,(self.root/'job.json').read_bytes())

    def test_queue_outage_rejects_submission(self):
        with patch.object(pc.subprocess,'run',return_value=subprocess.CompletedProcess([],1,'','trqauthd')):
            with self.assertRaisesRegex(ValueError,'调度服务无法核实'):pc.queue_guard(self.root)

    def test_active_recorded_job_rejected(self):
        (self.root/'提交记录.json').write_text(json.dumps([{'job_id':'123.node1'}]))
        results=[subprocess.CompletedProcess([],0,'user queue',''),subprocess.CompletedProcess([],0,'job_state = R','')]
        with patch.object(pc.subprocess,'run',side_effect=results):
            with self.assertRaisesRegex(ValueError,'已有作业活动中'):pc.queue_guard(self.root)

    def test_view_does_not_change_original_files(self):
        stage=self.root/'raw/A/conformers/c0001/opt';stage.mkdir(parents=True)
        inp=stage/'job.gjf';log=stage/'job.log';inp.write_text('# opt\n');log.write_text('incomplete output\n')
        originals={p:hashlib.sha256(p.read_bytes()).hexdigest() for p in [inp,log,self.root/'inputs/A.pdb']}
        z=pv.collect(self.root,True)
        self.assertGreaterEqual(z['files'],3)
        self.assertEqual((self.root/'过程文件/01_Gaussian/A/c0001/01_优化/job.log').read_text(),'incomplete output\n')
        self.assertTrue((self.root/'计算总览.html').exists())
        self.assertEqual(originals,{p:hashlib.sha256(p.read_bytes()).hexdigest() for p in originals})


if __name__=='__main__':unittest.main()
