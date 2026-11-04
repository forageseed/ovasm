import importlib.util
from pathlib import Path
import unittest

P = Path(__file__).resolve().parents[1] / "scripts/assembly_benchmarks/experimental_branch_posterior.py"
spec = importlib.util.spec_from_file_location("branch_posterior", P)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def dataset():
    return {"candidates": [{"pairings": [
        {"path": "u1+ r+ u2+", "read_counts": {"target": 9, "panel_a": 8, "panel_b": 2}},
        {"path": "u1+ r+ u3+", "read_counts": {"target": 1, "panel_a": 2, "panel_b": 8}},
    ]}]}


class ConditionalModelTest(unittest.TestCase):
    def test_empirical_prior_and_capped_target_update(self):
        out = module.compute(dataset(), "target", 2, 5, 1)
        self.assertEqual(len(out["groups"]), 1)
        g = out["groups"][0]
        self.assertEqual(g["panel_runs_used"], ["panel_a", "panel_b"])
        self.assertAlmostEqual(g["target_total_effective"], 5)
        a, b = g["alternatives"]
        self.assertAlmostEqual(a["prior_alpha"], 1.5)
        self.assertAlmostEqual(b["prior_alpha"], 1.5)
        self.assertAlmostEqual(a["conditional_mean"], 0.75)
        self.assertAlmostEqual(b["conditional_mean"], 0.25)
        self.assertAlmostEqual(a["conditional_variance"], .75*.25/9)

    def test_uncapped_iid_and_insufficient_panel(self):
        out = module.compute(dataset(), "target", 2, 0, 100)
        g = out["groups"][0]
        self.assertEqual(g["prior_status"], "insufficient_panel_runs_base_prior_only")
        self.assertAlmostEqual(g["alternatives"][0]["conditional_mean"], 9.5/11)
        self.assertEqual(g["applied_prior_strength"], 0)

    def test_zero_strength_does_not_claim_panel_was_used(self):
        g = module.compute(dataset(), "target", 0, 5, 1)["groups"][0]
        self.assertEqual(g["prior_status"], "panel_mean_computed_zero_strength")
        self.assertEqual(g["panel_runs_used"], [])
        self.assertEqual(g["panel_runs_eligible"], ["panel_a", "panel_b"])

    def test_duplicate_pairing_disagreement_rejected(self):
        d = dataset()
        d["candidates"].append({"pairings": [{"path": "u1+ r+ u2+", "read_counts": {"target": 1}}]})
        with self.assertRaises(ValueError):
            module.compute(d, "target", 2, 0, 1)

    def test_missing_is_not_treated_as_zero(self):
        d = dataset()
        del d["candidates"][0]["pairings"][1]["read_counts"]["target"]
        with self.assertRaises(ValueError):
            module.compute(d, "target", 2, 0, 1)
        d = dataset()
        del d["candidates"][0]["pairings"][1]["read_counts"]["panel_b"]
        g = module.compute(d, "target", 2, 0, 1)["groups"][0]
        self.assertEqual(g["panel_runs_eligible"], ["panel_a"])
        self.assertEqual(g["panel_runs_used"], [])
        self.assertEqual(g["prior_status"], "insufficient_panel_runs_base_prior_only")

    def test_reverse_complement_duplicate_is_rejected(self):
        d = dataset()
        d["candidates"].append({"pairings": [{"path": "u2- r- u1-", "read_counts": {"target": 9, "panel_a": 8, "panel_b": 2}}]})
        with self.assertRaises(ValueError):
            module.compute(d, "target", 2, 0, 1)

    def test_invalid_target_and_nonfinite_strength_rejected(self):
        with self.assertRaises(ValueError):
            module.compute(dataset(), "missing", 2, 0, 1)
        with self.assertRaises(ValueError):
            module.compute(dataset(), "target", float("nan"), 0, 1)


if __name__ == "__main__":
    unittest.main()
