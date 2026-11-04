import importlib.util
from pathlib import Path
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/assembly_benchmarks/compare_competing_alignments.py"
spec = importlib.util.spec_from_file_location("compare_competing_alignments", SCRIPT)
model = importlib.util.module_from_spec(spec)
spec.loader.exec_module(model)

def paf(name, target, score, start=0, end=150):
    n=end-start
    return "\t".join(map(str,[name,150,start,end,"+",target,1000,0,n,n,n,60,f"AS:i:{score}"]))+"\n"

class CompetingAlignments(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.path=Path(self.tmp.name)/"read.paf"
    def tearDown(self):
        self.tmp.cleanup()
    def test_nuclear_ambiguous_and_missing_are_not_origin_calls(self):
        self.path.write_text(paf("r1","mito",280)+paf("r1","nuclear",300)+paf("r2","mito",290)+paf("r2","nuclear",300)+paf("r3","mito",300))
        result=model.compare(self.path,["mito"],["nuclear"])
        self.assertEqual(result["summary"],{"nuclear_score_higher":1,"score_ambiguous":1,"missing_qualifying_competitor":1})
    def test_different_alignment_lengths_not_directly_compared(self):
        self.path.write_text(paf("r1","mito",300)+paf("r1","nuclear",280,0,140))
        result=model.compare(self.path,["mito"],["nuclear"])
        self.assertEqual(result["reads"][0]["status"],"alignment_span_incomparable")
    def test_ambiguous_target_sets_and_missing_scores_rejected(self):
        with self.assertRaisesRegex(ValueError,"disjoint"):
            model.compare(self.path,["same"],["same"])
        self.path.write_text(paf("r1","mito",280).replace("AS:i:280", "NM:i:2"))
        with self.assertRaisesRegex(ValueError,"AS:i"):
            model.compare(self.path,["mito"],["nuclear"])

if __name__=="__main__":
    unittest.main()
