import importlib.util
from pathlib import Path
import subprocess
import tempfile
import json
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
    def test_failed_summary_and_command_are_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp); root=p/"runs"; refs=p/"refs"
            sample=root/"default"/"Lathyrus_linifolius_PRJEB76347"
            (sample/"ovasm").mkdir(parents=True)
            (refs/"drLatLini1").mkdir(parents=True)
            (sample/"ovasm/summary.json").write_text(json.dumps({"failed":["plastid"]}))
            command='Command being timed: "ovasm run hifi --reads original.fastq.gz"\nExit status: 1\n'
            (sample/"time.txt").write_text(command)
            output=p/"report.json"
            with patch.object(m.sys,"argv",["compare",str(root),str(refs),"minimap2",str(output)]):m.main()
            row=json.loads(output.read_text())["samples"][0]
            self.assertEqual(row["status"],"partial")
            self.assertEqual(row["failed"],["plastid"])
            self.assertEqual(row["command_record_text"],command)
if __name__ == "__main__":unittest.main()
