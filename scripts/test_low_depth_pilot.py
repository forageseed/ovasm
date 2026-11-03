import importlib.util
from pathlib import Path
import os
import tempfile
import unittest
spec=importlib.util.spec_from_file_location("pilot",Path(__file__).with_name("run_low_depth_pilot.py"))
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class RunnerTests(unittest.TestCase):
    def config(self):return dict(cpus=sorted(os.sched_getaffinity(0))[:1],memory_bytes=64*1024**3,timeout_seconds=.1,reserve_bytes=0,oatk="/bin/true")
    def test_nonzero_retained(self):
        with tempfile.TemporaryDirectory() as tmp:
            result=m.execute(["/bin/false"],Path(tmp),self.config(),False)
            self.assertEqual(result["exit_code"],1)
    def test_timeout_kills_group(self):
        with tempfile.TemporaryDirectory() as tmp:
            result=m.execute(["/bin/sleep","10"],Path(tmp),self.config(),False)
            self.assertEqual(result["termination"],"timeout")
            self.assertLess(result["wall_seconds"],5)
if __name__=="__main__":unittest.main()
