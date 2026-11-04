import copy
import importlib.util
from pathlib import Path
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/assembly_benchmarks/compare_pan_branch_models.py"
spec = importlib.util.spec_from_file_location("compare_pan_branch_models", SCRIPT)
model = importlib.util.module_from_spec(spec)
spec.loader.exec_module(model)

def example():
    return {"candidates": [{"pairings": [
        {"path": "a+ r+ b+", "read_counts": {"s1": 20, "s2": 80, "s3": 10}},
        {"path": "a+ r+ c+", "read_counts": {"s1": 80, "s2": 20, "s3": 90}},
    ]}]}

class ComparePanBranchModels(unittest.TestCase):
    def test_heldout_counts_do_not_change_its_predictions(self):
        data = example()
        first = model.evaluate(data, 10)
        altered = copy.deepcopy(data)
        altered["candidates"][0]["pairings"][0]["read_counts"]["s1"] = 95
        altered["candidates"][0]["pairings"][1]["read_counts"]["s1"] = 5
        second = model.evaluate(altered, 10)
        a = next(e for e in first["events"] if e["heldout"] == "s1")
        b = next(e for e in second["events"] if e["heldout"] == "s1")
        for name in model.MODELS:
            self.assertEqual(a["metrics"][name]["predicted"], b["metrics"][name]["predicted"])
        self.assertNotEqual(a["metrics"]["sample_mean_fraction"]["cross_entropy"], b["metrics"]["sample_mean_fraction"]["cross_entropy"])

    def test_brier_is_one_hot_outcome_score_not_fraction_error(self):
        result = model.evaluate(example(), 10)
        event = next(e for e in result["events"] if e["heldout"] == "s1")
        uniform = event["metrics"]["uniform"]
        self.assertAlmostEqual(uniform["brier"], 0.5)
        self.assertAlmostEqual(uniform["squared_fraction_error"], 0.18)

    def test_reverse_orientation_is_one_pairing_and_inconsistency_rejected(self):
        data = example()
        duplicate = {"path": "b- r- a-", "read_counts": {"s1": 20, "s2": 80, "s3": 10}}
        data["candidates"][0]["pairings"].append(duplicate)
        self.assertEqual(model.evaluate(data, 10)["event_count"], 3)
        duplicate["read_counts"]["s1"] = 21
        with self.assertRaisesRegex(ValueError, "inconsistent counts"):
            model.evaluate(data, 10)

    def test_pseudocount_is_checked_and_changes_predictions(self):
        with self.assertRaisesRegex(ValueError, "pseudocount"):
            model.evaluate(example(), 10, 0)
        a = model.evaluate(example(), 10, 0.5)
        b = model.evaluate(example(), 10, 10.0)
        e1 = next(e for e in a["events"] if e["heldout"] == "s1")
        e2 = next(e for e in b["events"] if e["heldout"] == "s1")
        self.assertNotEqual(e1["metrics"]["sample_mean_fraction"]["predicted"], e2["metrics"]["sample_mean_fraction"]["predicted"])

    def test_missing_and_negative_counts_are_rejected(self):
        data = example()
        del data["candidates"][0]["pairings"][0]["read_counts"]["s2"]
        with self.assertRaisesRegex(ValueError, "same three or more runs"):
            model.evaluate(data, 10)
        data = example()
        data["candidates"][0]["pairings"][0]["read_counts"]["s2"] = -1
        with self.assertRaisesRegex(ValueError, "nonnegative integers"):
            model.evaluate(data, 10)

if __name__ == "__main__":
    unittest.main()
