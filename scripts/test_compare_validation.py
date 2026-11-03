import importlib.util
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch
spec=importlib.util.spec_from_file_location("comparison",Path(__file__).with_name("compare_validation.py"))
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class ComparisonTests(unittest.TestCase):
    def test_aligner_failure_is_not_zero_coverage(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp); f=p/"a.fa"; f.write_text(">a\nACGT\n")
            with patch.object(m.subprocess,"run",side_effect=subprocess.CalledProcessError(1,"minimap2")):
                with self.assertRaises(subprocess.CalledProcessError):m.compare(f,f,"minimap2",p)
    def test_empty_root_is_explicit_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            with patch.object(m.sys,"argv",["compare",str(p),str(p),"minimap2",str(p/"report.json")]):
                with self.assertRaisesRegex(ValueError,"no sample"):m.main()
if __name__ == "__main__":unittest.main()
