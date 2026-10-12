"""The status observer must count live receipts without a final conformer manifest."""
import hashlib
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import progress_status as ps


class LiveProgressChecks(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.cfg={'output':'raw','methods':['01','02','03'],'candidates':[{'id':'A'}]}
        self.write('job.json',self.cfg)

    def write(self,name,value):
        p=self.root/name;p.parent.mkdir(parents=True,exist_ok=True)
        p.write_text(json.dumps(value),encoding='utf-8')
        return p

    def receipt(self,cid,stage,status=None):
        row={'stage':stage}
        if status:row['status']=status
        return self.write('raw/A/conformers/'+cid+'/'+stage+'/receipt.json',row)

    def test_live_running_candidate_without_manifest(self):
        self.write('raw/A/state.json',{'status':'running'})
        d=self.root/'raw/A/sampling';d.mkdir()
        (d/'gfn2.xyz').write_text('2\n-1.0\nH 0 0 0\nH 0 0 1\n'*3)
        for c in ('c0001','c0002'):
            self.receipt(c,'opt1');self.receipt(c,'freq');self.receipt(c,'sp')
        self.receipt('c0001','nmr01');self.receipt('c0001','nmr02')
        event={'event':'started','time':'2026-10-12 20:50:59','cwd':str(self.root/'raw/A/conformers/c0001/nmr03'),'log':str(self.root/'raw/A/conformers/c0001/nmr03/job.log')}
        f=self.root/'运行日志/software_events.jsonl';f.parent.mkdir();f.write_text(json.dumps(event)+'\n')
        row=ps.snapshot(self.root)['candidates'][0]
        self.assertEqual((row['gfn2'],row['dft_inputs'],row['optimized'],row['positive_frequencies'],row['single_points']),(3,2,2,2,2))
        self.assertEqual(row['nmr_tasks'],2)
        self.assertIsNone(row['nmr_target_tasks'])
        self.assertIn('c0001',row['current'])
        self.assertIn('NMR03',row['current'])
        self.assertIn('待汇总',ps.format_snapshot(ps.snapshot(self.root)))

    def test_rejected_frequency_and_opt_retries_not_double_counted(self):
        self.receipt('c0001','opt1');self.receipt('c0001','opt2')
        self.receipt('c0001','freq','rejected_frequency')
        row=ps.snapshot(self.root)['candidates'][0]
        self.assertEqual(row['optimized'],1)
        self.assertEqual(row['positive_frequencies'],0)
        self.assertEqual(row['rejected_frequencies'],1)

    def test_final_manifest_supplies_real_targets_and_coverage(self):
        self.write('raw/A/state.json',{'status':'complete'})
        self.write('raw/A/conformers.json',{'ensemble_total':22,'dft_selected':6,'nmr_ids':['c0001','c0003'],'coverage_within_computed_ensemble':0.9914})
        for c in ('c0001','c0003'):
            for m in ('01','02','03'):self.receipt(c,'nmr'+m)
        row=ps.snapshot(self.root)['candidates'][0]
        self.assertEqual((row['nmr_tasks'],row['nmr_target_tasks'],row['nmr_complete_conformers']),(6,6,2))
        self.assertAlmostEqual(row['coverage_percent'],99.14)
        self.assertIn('99.14%',ps.format_snapshot(ps.snapshot(self.root)))

    def test_observer_does_not_write_any_file(self):
        self.write('raw/A/state.json',{'status':'running'})
        self.receipt('c0001','freq')
        before={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in self.root.rglob('*') if p.is_file()}
        ps.snapshot(self.root)
        after={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in self.root.rglob('*') if p.is_file()}
        self.assertEqual(before,after)

    def test_partial_xyz_not_presented_as_complete_count(self):
        d=self.root/'raw/A/sampling';d.mkdir(parents=True)
        (d/'gfn2.xyz').write_text('2\n-1\nH 0 0 0\n')
        row=ps.snapshot(self.root)['candidates'][0]
        self.assertIsNone(row['gfn2'])
        self.assertTrue(row['warnings'])

    def test_unstarted_and_unreadable_records_are_distinct(self):
        self.assertEqual(ps.snapshot(self.root)['candidates'][0]['status'],'尚未开始')
        p=self.root/'raw/A/state.json';p.parent.mkdir(parents=True);p.write_text('{')
        row=ps.snapshot(self.root)['candidates'][0]
        self.assertEqual(row['status'],'状态记录暂不可读')
        self.assertTrue(row['warnings'])

    def test_partial_last_event_keeps_previous_valid_event(self):
        self.write('raw/A/state.json',{'status':'running'})
        events=self.root/'运行日志/software_events.jsonl';events.parent.mkdir()
        event={'event':'started','time':'2026-10-12 20:50:59','cwd':str(self.root/'raw/A/conformers/c0002/nmr02')}
        events.write_text(json.dumps(event)+'\n{"event":')
        self.assertIn('NMR02',ps.snapshot(self.root)['candidates'][0]['current'])

    def test_original_status_entry_uses_live_records(self):
        import project_console as pc
        self.write('raw/A/state.json',{'status':'running'})
        self.receipt('c0001','nmr01')
        stream=io.StringIO()
        with patch.object(pc,'collect',return_value={'candidates':[['A','running','','','','','']]}),contextlib.redirect_stdout(stream):
            pc.status(self.root)
        self.assertIn('NMR01=1',stream.getvalue())
        self.assertIn('待最终汇总',stream.getvalue())


if __name__=='__main__':unittest.main()
