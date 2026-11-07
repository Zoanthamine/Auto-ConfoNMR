"""Ensure the aggregate includes scientific auxiliary results without quantum bodies."""
import json
from pathlib import Path
import tempfile
import unittest

from workflow_log import build,ensemble_table


class WorkflowLogChecks(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.out=self.root/'过程文件'
        self.cfg={'output':'过程文件','candidates':[{'id':'A'}],'max_conformers':6,'clustering':{},'shermo':{'T':298.15}}
        (self.root/'job.json').write_text(json.dumps(self.cfg),encoding='utf-8')

    def tearDown(self):self.tmp.cleanup()

    def put(self,path,text):
        p=self.out/path;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(text,encoding='utf-8');return p

    def test_real_auxiliary_bodies_energy_and_final_weights_are_searchable(self):
        md=self.put('A/sampling/md01.log','MD_RUN完整动力学')
        self.put('A/sampling/gfn2/isostat.log','CLUSTER_RAW_TABLE')
        self.put('A/sampling/gfn2/cluster.xyz','1\nEnergy = -10.0 a.u. #Cluster: 1\nH 0 0 0\n1\nEnergy = -9.999 a.u. #Cluster: 2\nH 0 0 0\n')
        self.put('A/shermo/c0001.log','THERMOCHEMISTRY_COMPLETE')
        self.put('module_exports/conformers/A/populations.csv','conformer,percent\nc0001,100\n')
        self.put('A/conformers.json',json.dumps({'records':[{'id':'c0001','gibbs_hartree':-10,'weight':1.0}]}))
        original=md.read_bytes();result=build(self.root)
        text=Path(result['path']).read_text(encoding='utf-8')
        for marker in ['MD_RUN完整动力学','CLUSTER_RAW_TABLE','THERMOCHEMISTRY_COMPLETE','c0001,100','最终自由能比例','0.627509']:
            self.assertIn(marker,text)
        self.assertEqual(md.read_bytes(),original)
        build(self.root)
        text=Path(result['path']).read_text(encoding='utf-8')
        self.assertEqual(text.count('MD_RUN完整动力学'),1)

    def test_quantum_body_excluded_but_helpers_kept_and_old_names_recognized(self):
        self.put('A/sampling/resp_gas/Gau_RESP_gas_A.out','SECRET_GAUSSIAN_BODY')
        self.put('A/sampling/resp_gas/Gau_RESP_gas_A_Multiwfn.log','MULTIWFN_RESULTS')
        self.put('A/conformers/c0001/opt1/job.log','LEGACY_GAUSSIAN_BODY')
        self.put('A/conformers/c0001/sp/Orca_sp_A_C0001.log','ORCA_BODY')
        text=Path(build(self.root)['path']).read_text(encoding='utf-8')
        self.assertIn('MULTIWFN_RESULTS',text)
        self.assertIn('Orca_sp_A_C0001.log',text)
        for marker in ['SECRET_GAUSSIAN_BODY','LEGACY_GAUSSIAN_BODY','ORCA_BODY']:self.assertNotIn(marker,text)

    def test_incomplete_ensemble_is_not_reported_as_complete(self):
        path=self.put('A/sampling/gfn2/cluster.xyz','2\n-10.0\nH 0 0 0\n')
        with self.assertRaises(ValueError):ensemble_table(path,6)
        result=build(self.root)
        self.assertTrue(result['warnings'])
        self.assertNotIn('构象数=1',Path(result['path']).read_text(encoding='utf-8'))

    def test_binary_temporary_and_static_resource_contents_not_merged(self):
        self.put('A/sampling/OPTIM/1/run.log','IN_PROGRESS_XTB_TEMP')
        self.put('A/sampling/molecule.xtc','BINARY_TRAJECTORY')
        self.put('A/sampling/bonded_param.dat','STATIC_FORCE_FIELD')
        resource=self.root/'workflow/resources/bonded_param.dat';resource.parent.mkdir(parents=True);resource.write_text('STATIC_FORCE_FIELD')
        text=Path(build(self.root)['path']).read_text(encoding='utf-8')
        self.assertIn('molecule.xtc',text)
        for marker in ['IN_PROGRESS_XTB_TEMP','BINARY_TRAJECTORY','STATIC_FORCE_FIELD']:self.assertNotIn(marker,text)


if __name__=='__main__':unittest.main()
