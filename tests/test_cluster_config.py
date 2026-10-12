"""Portability and resource checks without launching real chemistry programs."""
import copy
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import cluster_config as cc
import nmrflow as nf
import project_console as pc


class ClusterConfigurationChecks(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.base = Path(self.tmp.name)
        self.site = json.loads((cc.ROOT/'config/cluster.example.json').read_text())
        for key in ('projects_root','scratch_root','resource_bundle'):
            self.site[key] = str(self.base/key)
        self.site['python'] = sys.executable
        self.site['software'] = {k:str(self.base/'software'/k) for k in cc.SOFTWARE}
        self.site['gaussian_profile'] = ''
        self.site['scheduler']['queue'] = 'test_queue'
        self.path = self.base/'site.json'

    def tearDown(self):
        self.tmp.cleanup()

    def load(self):
        self.path.write_text(json.dumps(self.site), encoding='utf-8')
        return cc.load_site(self.path)

    def resources(self):
        src = Path(self.site['resource_bundle'])
        for name in list(cc.EXTERNAL_RESOURCES)+['upstream/'+n for n in cc.EXTERNAL_WORKBOOKS]:
            path = src/name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(('fixture '+name).encode())

    def test_missing_site_never_falls_back_to_author_paths(self):
        with self.assertRaisesRegex(ValueError, '缺少集群配置'):
            cc.load_site(self.base/'missing.json')

    def test_placeholders_and_unknown_fields_fail(self):
        self.site['scheduler']['queue'] = 'CHANGE_ME'
        with self.assertRaises(ValueError): self.load()
        self.site['scheduler']['queue'] = 'queue'
        self.site['unexpected'] = 3
        with self.assertRaisesRegex(ValueError, '字段'): self.load()

    def test_workbooks_required_and_installed_byte_exact(self):
        self.resources()
        self.assertEqual(cc.missing_resources(self.site), [])
        target = self.base/'workflow'
        cc.install_resources(self.site, target)
        for name in cc.EXTERNAL_WORKBOOKS:
            self.assertEqual((target/'upstream'/name).read_bytes(),
                             (Path(self.site['resource_bundle'])/'upstream'/name).read_bytes())
        (Path(self.site['resource_bundle'])/'upstream/STS.xlsx').unlink()
        with self.assertRaisesRegex(ValueError, 'STS.xlsx'): cc.install_resources(self.site, target)

    def test_portable_resource_paths_and_xtb_data(self):
        self.resources()
        source = Path(self.site['resource_bundle'])
        (source/'share/xtb').mkdir(parents=True)
        (source/'sobtop.ini').write_text(
            '  nthreads= 4\n  k_method= 2\n  Multiwfn_cmd= "old"\n  OpenBabel_cmd= "old"\n')
        target = self.base/'workflow'
        cc.install_resources(self.site, target)
        actual = (target/'resources/sobtop.ini').read_text()
        self.assertIn('nthreads= 4', actual)
        self.assertIn('k_method= 2', actual)
        self.assertIn('Multiwfn_cmd= "'+self.site['software']['multiwfn']+'"', actual)
        self.assertIn('OpenBabel_cmd= "'+self.site['software']['obabel']+'"', actual)
        self.assertNotIn('"old"', actual)
        self.assertIn('export XTBPATH=', cc.runtime_text(self.site))

    def test_scheduler_cores_memory_coupled(self):
        site = self.load()
        site['scheduler']['node_resource'] = 'select=1:ncpus={cores}:mem={memory_gb}gb'
        site['scheduler']['extra_qsub_args'] = ['-A','account']
        args = cc.scheduler_arguments(site, {'cores':8,'memory_gb':24})
        self.assertEqual(args, ['-q','test_queue','-l','select=1:ncpus=8:mem=24gb',
                                '-l','walltime=72:00:00','-A','account'])

    def test_runtime_quotes_paths_and_modules(self):
        site = self.load()
        site['python'] = '/migration/env/bin/python'
        site['software']['multiwfn'] = '/migration/Multiwfn_noGUI'
        site['modules'] = ['gromacs/2024']
        site['path_prepend'] = ["/a path/has'quote"]
        text = cc.runtime_text(site)
        self.assertIn('module load gromacs/2024', text)
        self.assertIn("'\"'\"'", text)
        self.assertNotIn('huangjt', text)
        site['modules'] = []
        self.assertNotIn('module load', cc.runtime_text(site))

    def test_prepare_missing_workbook_does_not_create_project(self):
        source = self.base/'A.pdb';source.write_text('structure')
        with patch.object(pc, 'load_site', return_value=self.site):
            with self.assertRaisesRegex(ValueError, '外部资源'):
                pc.prepare(self.base/'new_project', [source])
        self.assertFalse((self.base/'new_project').exists())

    def test_wrong_allocation_fails_before_any_program(self):
        cfg = json.loads((nf.ROOT/'job.example.json').read_text())
        cfg['cores'] = 8
        nodefile = self.base/'nodes'
        for lines in ('node1\nnode2\n','node1\n'*4):
            nodefile.write_text(lines)
            with patch.dict(os.environ, {'PBS_JOBID':'test','PBS_NODEFILE':str(nodefile)}), \
                    patch.object(nf, 'run_candidate') as run:
                with self.assertRaisesRegex(ValueError, 'CPU slots'):
                    nf.execute('run',cfg,self.base/'out',self.base)
                run.assert_not_called()


if __name__ == '__main__':
    unittest.main()
