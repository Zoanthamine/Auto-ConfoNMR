"""Check selected science reaches real inputs, output labels and restart guards."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import nmrflow as nf
import project_console as pc
from protocol_options import apply_choices, cluster_options, method_map, nmr_text, template_method


class ProtocolChecks(unittest.TestCase):
    def cfg(self):return json.loads((nf.ROOT/'job.example.json').read_text(encoding='utf-8'))

    def test_new_order_and_calibration_mismatch_rejected(self):
        cfg=self.cfg();nf.validate_config(cfg)
        self.assertEqual([method_map(cfg)[m][2] for m in cfg['methods']],['dp4plus','dp4plus','sts'])
        incomplete=json.loads(json.dumps(cfg));incomplete.pop('nmr_settings')
        with self.assertRaisesRegex(ValueError,'缺少显式'):nf.validate_config(incomplete)
        cfg['postprocessing'][1]['method']='03'
        with self.assertRaisesRegex(ValueError,'parameter does not match'):nf.validate_config(cfg)

    def test_profiles_change_actual_gaussian_input_and_stage_hash(self):
        geom=[['C',0,0,0],['H',1,0,0],['H',0,1,0],['H',0,0,1],['H',-1,0,0]]
        from test_nmrflow import gaussian
        for profile,model,keyword in [('legacy_pcm_g09','IEFPCM',True),('pcm_g16','IEFPCM',False),('smd_g16','SMD',False)]:
            cfg=self.cfg();apply_choices(cfg,profile)
            for mid,functional in [('01','b3lyp/6-31g(d)'),('02','mpw1pw91/6-311+g(d,p)'),('03','wb97xd/6-31g(d)')]:
                with self.subTest(profile=profile,mid=mid),tempfile.TemporaryDirectory() as d:
                    def fake(argv,cwd,log,env=None,stdin=None):
                        self.assertIn(functional,stdin.lower());self.assertIn('SCRF=('+model+',Solvent=Chloroform)',stdin)
                        self.assertEqual('G09Defaults' in stdin,keyword)
                        Path(log).write_text(gaussian(geom))
                    with patch.object(nf,'command',fake):nf.quantum(nf.ROOT/'resources'/('06_tem'+mid+'_nmr.gjf'),geom,cfg,Path(d),'nmr')
                    changed=json.loads(json.dumps(cfg));apply_choices(changed,'smd_g16' if profile!='smd_g16' else 'legacy_pcm_g09')
                    with self.assertRaisesRegex(ValueError,'Changed input/output'):
                        nf.quantum(nf.ROOT/'resources'/('06_tem'+mid+'_nmr.gjf'),geom,changed,Path(d),'nmr')

    def test_clustering_commands_explicitly_use_all_stage_values(self):
        cfg=self.cfg()
        for name,expected in [('coarse',[(.5,.5,300),(.3,.2,300),(.3,.2,300)]),('fine',[(.2,.2,350),(.2,.1,350),(.2,.1,350)])]:
            apply_choices(cfg,clustering_profile=name)
            for stage,values in zip(['GFN1','GFN2','DFT'],expected):
                argv=cluster_options(cfg,stage);actual=dict(zip(argv[::2],argv[1::2]))
                self.assertEqual(tuple(actual[k] for k in ['-Edis','-Gdis','-T']),values)
                self.assertEqual(actual['-Nout'],100);self.assertEqual(actual['-Eout'],100)

    def test_historic_mapping_and_smd_routes(self):
        cfg=self.cfg()
        for k in ['method_layout','nmr_settings','nmr_profile','clustering_profile']:cfg.pop(k,None)
        cfg['clustering']={'Edis':.2,'GFN1_Gdis':.2,'Gdis':.1}
        self.assertEqual(method_map(cfg)['02'][0],'wB97XD');self.assertEqual(template_method('02',cfg),'03')
        text=(nf.ROOT/'resources/06_tem03_nmr.gjf').read_text()
        rendered=nmr_text(text,cfg)
        self.assertIn('SCRF=(SMD,Solvent=Chloroform)',rendered)
        self.assertNotIn('G09Defaults',rendered)
        historic='# NMR wB97XD/6-31G(d) scrf(SMD,solvent=chloroform)\n'
        self.assertEqual(nmr_text(historic,cfg),historic)
        self.assertEqual(cluster_options(cfg,'GFN2'),['-Edis',.2,'-Gdis',.1,'-T',350])

    def test_ini_choices_change_before_run_but_block_after_run(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);cfg=self.cfg();cfg.update(output='raw',candidates=[{'id':'A','input':'inputs/A.pdb'}])
            (root/'inputs').mkdir();(root/'inputs/A.pdb').write_text('COMPND A\n')
            (root/'job.json').write_text(json.dumps(cfg));ini=root/'项目配置.ini';ini.write_text(pc.settings_text(cfg),encoding='utf-8')
            ini.write_text(ini.read_text(encoding='utf-8').replace('clustering_profile = coarse','clustering_profile = fine'),encoding='utf-8')
            _,cfg,_=pc.sync(root);self.assertEqual(cfg['clustering']['GFN2']['Gdis'],.1)
            self.assertIn('mPW1PW91', (root/'参数预览/NMR_02_模板_含坐标占位符.gjf').read_text(encoding='utf-8'))
            resolved=json.loads(json.dumps(cfg));resolved['output']=str(root/'raw');resolved['candidates'][0]['input']=str(root/'inputs/A.pdb')
            state=root/'raw/A/state.json';state.parent.mkdir(parents=True);state.write_text(json.dumps({'identity':{'config':nf.calculation_config(resolved)}}))
            old=(root/'job.json').read_bytes()
            ini.write_text(ini.read_text(encoding='utf-8').replace('nmr_profile = legacy_pcm_g09','nmr_profile = pcm_g16'),encoding='utf-8')
            with self.assertRaisesRegex(ValueError,'计算配置发生变化'):pc.sync(root)
            self.assertEqual(old,(root/'job.json').read_bytes())

    def test_project_cannot_be_created_inside_workflow(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);program=root/'program';program.mkdir();project=program/'nested_project'
            with patch.object(pc,'ROOT',program),self.assertRaisesRegex(ValueError,'程序目录外'):
                pc.prepare(project,[])
            self.assertFalse(project.exists())


if __name__=='__main__':unittest.main()
