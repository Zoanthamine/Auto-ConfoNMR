"""V2 numerical, routing, restart and data-integrity checks."""
import configparser
import contextlib
import io
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import nmrflow as nf
import module_options as mo
import conformer_exports as ce
import coupling as jc
import ecd
import project_console as pc


def td_text():
    return '''Excited State   1: Singlet-A  4.0000 eV 309.96 nm f=0.0100
Excited State   2: Singlet-A  5.0000 eV 247.97 nm f=0.0200
 Rotatory Strengths (R) in cgs (10**-40 erg-esu-cm/Gauss)
       state          XX          YY          ZZ    R(velocity)    E-M Angle
         1         1.0000      0.0000      0.0000      1.0000       90.00
         2        -2.0000      0.0000      0.0000     -2.0000       90.00
 END
 Rotatory Strengths (R) in cgs (10**-40 erg-esu-cm/Gauss)
       state          XX          YY          ZZ    R(length)
         1         1.1000      0.0000      0.0000      1.1000
         2        -2.2000      0.0000      0.0000     -2.2000
 END
'''


class ModuleChecks(unittest.TestCase):
    def cfg(self):return json.loads((nf.ROOT/'job.example.json').read_text())

    def test_nmr_without_experiment_refreshes_raw_exports(self):
        cfg=self.cfg();cfg['assignments']='assignments.csv'
        p=configparser.ConfigParser()
        p['对比']={key:'' for key in ('experiment','assignments','upper_experiment','upper_assignments','lower_experiment','lower_assignments')}
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);job=root/'job.json'
            with patch.object(pc,'sync',return_value=(job,cfg,p)), \
                 patch.object(pc,'load_config',return_value=(job,cfg,root/'过程文件')), \
                 patch.object(pc.subprocess,'run') as run, \
                 patch.object(pc,'queue_guard'), \
                 patch.object(pc,'collect') as collect, contextlib.redirect_stdout(io.StringIO()):
                pc.compare(root)
            run.assert_called_once_with([pc.sys.executable,str(pc.ROOT/'nmrflow.py'),'post',str(job)],check=True)
            collect.assert_called_once_with(root,True)

    def test_defaults_keep_nmr_and_no_extra_quantum(self):
        cfg=self.cfg();self.assertEqual(mo.options(cfg)['base'],'nmr')
        self.assertFalse(mo.options(cfg)['coupling']);self.assertFalse(mo.options(cfg)['ecd']);mo.validate(cfg)

    def test_quick_incompatible_with_properties_rejected(self):
        cfg=self.cfg();cfg['modules']=dict(mo.DEFAULT_MODULES,base='conformers',conformer_weight='gfn2_energy',ecd=True)
        with self.assertRaisesRegex(ValueError,'standalone'):mo.validate(cfg)

    def test_analysis_edit_is_not_a_new_calculation(self):
        cfg=self.cfg();original=nf.calculation_config(cfg)
        cfg['module_analysis']={'top_n':12,'ecd_shift_nm':-3.0}
        self.assertEqual(original,nf.calculation_config(cfg))
        cfg['modules']=dict(mo.DEFAULT_MODULES,ecd=True)
        self.assertNotEqual(original,nf.calculation_config(cfg))

    def test_fast_population_coordinate_count_and_coverage(self):
        geom=[['H',0,0,0],['H',0,0,.74]]
        blocks=[('-1.0',geom),('-0.999',geom),('-0.998',geom)]
        rows=ce.electronic_records(blocks,298.15)
        expected=np.exp(-np.array([0,.001,.002])*627.509474/(.00198720425864083*298.15));expected/=sum(expected)
        np.testing.assert_allclose([r['weight'] for r in rows],expected)
        manifest={'candidate':'A','records':rows,'weight_model':'gfn2_electronic_energy','config':self.cfg()}
        with tempfile.TemporaryDirectory() as d:
            result=ce.export_candidate(manifest,Path(d),2)
            self.assertEqual(len(nf.xyz_read(Path(d)/'top_conformers.xyz')),2)
            self.assertAlmostEqual(result['selected_coverage'],sum(expected[:2]));self.assertLess(result['selected_coverage'],1)

    def test_fast_xyz_requires_declared_model_and_skips_all_quantum(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);source=root/'input.xyz';source.write_text(nf.xyz_text([['H',0,0,0],['H',0,0,.74]],'-1.0'))
            cfg=self.cfg();cfg['modules']=dict(mo.DEFAULT_MODULES,base='conformers',conformer_weight='gfn2_energy')
            candidate={'id':'A','input':str(source)}
            with patch.object(nf,'quantum',side_effect=AssertionError('Quantum called')):
                with self.assertRaisesRegex(ValueError,'input_energy_model'):nf.run_candidate(candidate,cfg,root/'raw/A')
                cfg['input_energy_model']='gfn2'
                nf.run_candidate(candidate,cfg,root/'raw2/A')
            self.assertEqual(json.loads((root/'raw2/A/conformers.json').read_text())['dft_selected'],0)

    def test_fc_block_columns_10_and_larger_keep_signs(self):
        n=12;lines=['Fermi Contact (FC) contribution to J (Hz):']
        expected={}
        for start in range(1,n+1,5):
            cols=list(range(start,min(n+1,start+5)));lines.append(' '.join(map(str,cols)))
            for i in range(start,n+1):
                values=[]
                for j in cols:
                    if j>i:continue
                    v=(-1 if (i+j)%2 else 1)*(i+j/100);expected[jc.pair_key(i,j)]=v
                    values.append('{:.6E}'.format(v).replace('E','D'))
                lines.append(str(i)+' '+' '.join(values))
        lines.append('Total nuclear spin-spin coupling J (Hz):')
        self.assertEqual(jc.parse_fc('\n'.join(lines),n),expected)
        with self.assertRaisesRegex(ValueError,'exactly one'):jc.parse_fc('\n'.join(lines*2),n)
        with self.assertRaisesRegex(ValueError,'Incomplete'):jc.parse_fc('Fermi Contact (FC) contribution to J (Hz):\n1\n1 0.0\nEND',2)

    def test_methyl_signed_average_before_magnitude(self):
        matrix={'1,4':-3.0,'2,4':1.0,'3,4':1.0}
        self.assertAlmostEqual(jc.group_average(matrix,[1,2,3],[4]),-1/3)
        with self.assertRaisesRegex(ValueError,'overlapping'):jc.group_average(matrix,[1,2],[2,3])

    def test_vicinal_graph_distance(self):
        graph={1:{2},2:{1,3},3:{2,4},4:{3}}
        self.assertEqual(jc.distance(graph,1,4),3);self.assertEqual(jc.distance(graph,1,3),2)

    def test_dj_matches_independently_verified_scipy_reference(self):
        data={}
        for cid,offset in [('A',0),('B',3)]:
            data[cid]=[{'type':'C','signal':str(i),'experimental':x,'calculated':jc.PARAMETERS['TMS']['C']-x} for i,x in enumerate([20,60,150])]
            data[cid].append({'type':'J','signal':'J1','experimental':5.0,'calculated':jc.PARAMETERS['j_slope']*5+jc.PARAMETERS['j_intercept']+offset})
        result,details=jc.probabilities(data)
        # Independently evaluated with scipy.stats.t.sf((3/0.9509)/0.992,3).
        # Keeping the checked fixture avoids a new SciPy runtime dependency.
        tails=[.5,.025040303190157313]
        self.assertAlmostEqual(result['J']['A'],100*tails[0]/sum(tails),10)
        self.assertAlmostEqual(sum(result['all'].values()),100)

    def test_ecd_parser_both_gauges_and_duplicate_rejection(self):
        rows=ecd.parse_transitions(td_text());self.assertEqual(rows[1]['r_velocity'],-2)
        self.assertEqual(rows[1]['r_length'],-2.2)
        with self.assertRaisesRegex(ValueError,'Multiple'):ecd.parse_transitions(td_text()*2)
        with self.assertRaisesRegex(ValueError,'Incomplete'):ecd.parse_transitions(td_text().replace('         2        -2.2000      0.0000      0.0000     -2.2000\n',''))

    def test_ecd_formula_peak_fwhm_and_linear_weighting(self):
        width=.4;energy=4.0;rows=[{'energy_ev':energy,'r_velocity':10.0}]
        xs=ecd.HC_EV_NM/np.array([energy-width/2,energy,energy+width/2])
        curve=ecd.spectrum(rows,xs,width)
        self.assertAlmostEqual(curve[0]/curve[1],.5);self.assertAlmostEqual(curve[2]/curve[1],.5)
        sigma=width/(2*math.sqrt(math.log(2)))
        self.assertAlmostEqual(curve[1],energy*10/(22.96*sigma*math.sqrt(math.pi)))
        mirror=[{'energy_ev':energy,'r_velocity':-10.0}]
        np.testing.assert_allclose(ecd.spectrum(mirror,xs,width),-curve)

    def test_ecd_experiment_preserved_exactly(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'e.csv';p.write_text('wavelength_nm,delta_epsilon\n600,17\n550,-20\n200,30\n')
            before=p.read_bytes();values=ecd.experimental(p)
            self.assertEqual(values,[(600,17),(550,-20),(200,30)]);self.assertEqual(before,p.read_bytes())

    def test_ecd_does_not_inherit_nmr_pcm_g09(self):
        from test_nmrflow import gaussian
        cfg=self.cfg();cfg['modules']=dict(mo.DEFAULT_MODULES,ecd=True,ecd_states=2,ecd_solvent='acetonitrile')
        geom=[['H',0,0,0],['H',0,0,.74]]
        def fake(argv,cwd,log,env=None,stdin=None):
            self.assertIn('SMD,solvent=acetonitrile',stdin);self.assertIn('nstates=2',stdin)
            self.assertNotIn('G09Defaults',stdin);self.assertNotIn('IEFPCM',stdin)
            Path(log).write_text(gaussian(geom)+td_text())
        with tempfile.TemporaryDirectory() as d,patch.object(nf,'command',fake):
            nf.quantum(nf.ROOT/'resources/07_tem01_ecd.gjf',geom,cfg,Path(d),'ecd')

    def test_unfinished_candidates_cannot_postprocess(self):
        cfg=self.cfg()
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ValueError,'All candidates'):nf.execute('post',cfg,Path(d),Path(d))


if __name__=='__main__':unittest.main()
