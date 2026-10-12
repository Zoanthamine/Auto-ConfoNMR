"""Real files, identifiable names, quiet quantum log and checkpoint call chains."""
import contextlib
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import nmrflow as nf
import project_views as pv
import sampling as sp
from tests.test_nmrflow import gaussian


class DirectOutputChecks(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.cfg=json.loads((nf.ROOT/'job.example.json').read_text())
        self.cfg['output_layout']='direct_v1'
        self.geom=[['C',0,0,0],['H',1,0,0],['H',0,1,0],['H',0,0,1],['H',-1,0,0]]

    def tearDown(self):self.tmp.cleanup()

    def test_muted_quantum_keeps_full_bytes_and_progress_on_failure(self):
        log=self.root/'Gau_opt_001RR_C0001_A01.log';capture=io.StringIO()
        env=dict(os.environ,NMR_DETAIL_LOG='0')
        with patch.dict(os.environ,{'NMR_STREAM_LOG':'1','NMR_PROJECT_DIR':''}),contextlib.redirect_stdout(capture):
            with self.assertRaisesRegex(ValueError,'Command failed'):
                nf.command([sys.executable,'-c',"import sys;print('QUANTUM_DETAIL');print('ERR',file=sys.stderr);sys.exit(7)"],self.root,log,env=env)
        self.assertIn('QUANTUM_DETAIL',log.read_text());self.assertIn('ERR',log.read_text())
        self.assertNotIn('\nQUANTUM_DETAIL\n',capture.getvalue());self.assertNotIn('\nERR\n',capture.getvalue())
        self.assertIn('软件开始',capture.getvalue());self.assertIn('返回值 7',capture.getvalue())

    def test_quantum_checkpoint_names_and_resume(self):
        directory=self.root/'001RR/conformers/c0001/opt2'
        def fake(argv,cwd,log,env=None,stdin=None):
            self.assertIn('%chk=Gau_opt_001RR_C0001_A02.chk',stdin)
            self.assertEqual(env['NMR_DETAIL_LOG'],'0')
            Path(log).write_text(gaussian(self.geom,'opt'))
        with patch.object(nf,'command',side_effect=fake) as run:
            for _ in range(2):nf.quantum(nf.ROOT/'resources/03_tem2.gjf',self.geom,self.cfg,directory,'opt')
            self.assertEqual(run.call_count,1)
        self.assertTrue((directory/'Gau_opt_001RR_C0001_A02.gjf').is_file())
        self.assertFalse((directory/'job.log').exists())
        self.assertEqual(nf.artifact_stem(directory.parent/'nmr02',self.cfg,'Gau','nmr'),'Gau_nmr02_001RR_C0001')

    def test_orca_input_and_side_products_share_identifiable_stem(self):
        directory=self.root/'002SR/conformers/c0007/sp'
        def fake(argv,cwd,log,env=None,stdin=None):
            self.assertEqual(argv[1],'Orca_sp_002SR_C0007.inp')
            self.assertEqual(env['NMR_DETAIL_LOG'],'0')
            Path(log).write_text('FINAL SINGLE POINT ENERGY -40.1\nORCA TERMINATED NORMALLY\n')
            (Path(cwd)/(Path(argv[1]).stem+'.gbw')).write_bytes(b'wavefunction')
        with patch.object(nf,'command',side_effect=fake) as run:
            self.assertEqual(nf.orca_sp(self.geom,self.cfg,directory),-40.1)
            nf.orca_sp(self.geom,self.cfg,directory)
            self.assertEqual(run.call_count,1)
        self.assertTrue((directory/'Orca_sp_002SR_C0007.gbw').exists())

    def test_resp_formchk_multiwfn_file_chain_uses_new_stem(self):
        work=self.root/'001RR/sampling';work.mkdir(parents=True)
        def fake(argv,cwd,log,env=None,stdin=None):
            argv=[str(v) for v in argv];cwd=Path(cwd);Path(log).write_text('normal utility output')
            if argv[0]==self.cfg['software']['obabel']:
                (work/'input.xyz').write_text(nf.xyz_text(self.geom))
            elif argv[0]==self.cfg['software']['g16']:
                self.assertEqual(env['NMR_DETAIL_LOG'],'0')
                stem=Path(log).stem;self.assertIn('%chk='+stem+'.chk',stdin)
                self.assertIn('B3LYP/6-31g pop=MK IOp(6/33=2,6/42=6)',stdin)
                Path(log).write_text('Normal termination of Gaussian')
                (cwd/(stem+'.chk')).write_bytes(b'chk')
            elif argv[0]==self.cfg['software']['formchk']:
                self.assertTrue((cwd/argv[1]).is_file());(cwd/argv[2]).write_text('fchk')
            elif argv[0]==self.cfg['software']['multiwfn']:
                stem=Path(argv[1]).stem;self.assertIn(stem+'.out',stdin)
                (cwd/(stem+'.chg')).write_text('\n'.join(' '.join(map(str,a+[0])) for a in self.geom))
        with patch.object(sp,'command',side_effect=fake):sp.resp2(work/'molecule.mol2',self.cfg,work)
        self.assertTrue((work/'resp_solvent/Gau_RESP_solvent_001RR.fchk').exists())
        self.assertEqual(sum(json.loads((work/'resp2.json').read_text())['mixed']),0)

    def test_collect_real_large_files_no_symlinks_no_alias_copies(self):
        self.cfg.update(output='过程文件',candidates=[{'id':'001RR','input':'inputs/001RR.pdb'}])
        (self.root/'inputs').mkdir();(self.root/'inputs/001RR.pdb').write_text('input')
        (self.root/'job.json').write_text(json.dumps(self.cfg))
        stage=self.root/'过程文件/001RR/conformers/c0001/opt1';stage.mkdir(parents=True)
        log=stage/'Gau_opt_001RR_C0001_A01.log';log.write_text('still running')
        chk=stage/'Gau_opt_001RR_C0001_A01.chk'
        with chk.open('wb') as f:f.truncate(40*1024*1024)
        before=(chk.stat().st_ino,chk.stat().st_size)
        for copy in [False,True]:pv.collect(self.root,copy)
        self.assertEqual(before,(chk.stat().st_ino,chk.stat().st_size))
        self.assertFalse(any(p.is_symlink() for p in (self.root/'过程文件').rglob('*')))
        self.assertEqual(len(list((self.root/'过程文件').rglob('*.chk'))),1)
        self.assertIn('Gau_opt_001RR_C0001_A01.log',(self.root/'文件索引.csv').read_text(encoding='utf-8-sig'))


if __name__=='__main__':unittest.main()
