"""The same shared-directory gate must serialize kernel-capable and ENOLCK clients."""
import errno
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import nmrflow as nf


@unittest.skipIf(os.name=='nt','POSIX NFS kernel-lock regression')
class SharedLockChecks(unittest.TestCase):
    def test_enolck_still_has_cross_process_exclusion(self):
        import fcntl
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)
            with patch.object(fcntl,'flock',side_effect=OSError(errno.ENOLCK,'No locks available')):
                with nf.output_lock(out):
                    owner=json.loads((out/'.workflow.guard/owner.json').read_text())
                    self.assertEqual(owner['pid'],os.getpid())
                    self.assertEqual(owner['kernel_lock_errno'],errno.ENOLCK)
                    code="import sys;from pathlib import Path;sys.path.insert(0,"+repr(str(nf.ROOT))+");import nmrflow as nf\ntry:\n with nf.output_lock(Path("+repr(tmp)+")):print('UNSAFE_ENTRY')\nexcept ValueError:print('BUSY')"
                    r=subprocess.run([sys.executable,'-c',code],capture_output=True,text=True)
                    self.assertEqual(r.returncode,0,r.stderr)
                    self.assertEqual(r.stdout.strip(),'BUSY')
            self.assertFalse((out/'.workflow.guard').exists())

    def test_existing_guard_is_not_removed_or_stolen(self):
        import fcntl
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp);guard=out/'.workflow.guard';guard.mkdir()
            owner=guard/'owner.json';owner.write_text('{"token":"another-process"}')
            with patch.object(fcntl,'flock') as kernel:
                with self.assertRaisesRegex(ValueError,'Another workflow'):
                    with nf.output_lock(out):self.fail('Entered an existing guard')
                kernel.assert_not_called()
            self.assertEqual(owner.read_text(),'{"token":"another-process"}')

    def test_kernel_contention_cleans_only_own_gate(self):
        import fcntl
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)
            with patch.object(fcntl,'flock',side_effect=BlockingIOError(errno.EAGAIN,'busy')):
                with self.assertRaisesRegex(ValueError,'Another workflow'):
                    with nf.output_lock(out):self.fail('Entered a busy kernel lock')
            self.assertFalse((out/'.workflow.guard').exists())

    def test_unexpected_kernel_error_is_not_treated_as_contention(self):
        import fcntl
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)
            with patch.object(fcntl,'flock',side_effect=OSError(errno.EIO,'IO failed')):
                with self.assertRaisesRegex(ValueError,'lock service failed') as caught:
                    with nf.output_lock(out):self.fail('Ignored unknown kernel error')
            self.assertEqual(caught.exception.__cause__.errno,errno.EIO)
            self.assertFalse((out/'.workflow.guard').exists())


if __name__=='__main__':unittest.main()
