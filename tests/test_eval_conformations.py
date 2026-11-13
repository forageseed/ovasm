import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/eval_conformations.py"
MINIMAP2 = os.environ.get("MINIMAP2") or shutil.which("minimap2")


def rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


@unittest.skipUnless(MINIMAP2, "minimap2 not found (set MINIMAP2)")
class EvalConformationsTest(unittest.TestCase):
    """A 3 kb repeat X with entries a1, a2 and exits b1, b2. The molecules are a1-X-b1 and
    a2-X-b2 only, so two of the four traversals the graph allows are false. The molecules are
    circular through b1->a2 and b2->a1, which the graph lacks (missed links)."""

    def molecules(self, rng):
        s = {n: rand_seq(rng, 8000) for n in ("a1", "a2", "b1", "b2")}
        s["X"] = rand_seq(rng, 3000)
        # the exits share their first 900 bp (more than the 500 bp anchor), as in Geum rivale
        s["b2"] = s["b1"][:900] + s["b2"][900:]
        circle = s["a1"] + s["X"] + s["b1"] + s["a2"] + s["X"] + s["b2"]
        return s, circle

    def test_split_segments_are_compacted(self):
        """X stored as X1 (reverse strand) + X2 with a 50 bp overlap, a1 as two pieces."""
        rng = random.Random(7)
        s, circle = self.molecules(rng)
        rc = lambda x: x.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        seg = {"a1p": s["a1"][:5000], "a1q": s["a1"][5000:], "a2": s["a2"], "b1": s["b1"], "b2": s["b2"],
               "X1": rc(s["X"][:1550]), "X2": s["X"][1500:]}
        links = [("a1p", "+", "a1q", "+", 0), ("a1q", "+", "X1", "-", 0), ("a2", "+", "X1", "-", 0),
                 ("X1", "-", "X2", "+", 50), ("X2", "+", "b1", "+", 0), ("X2", "+", "b2", "+", 0)]
        gfa = "".join(f"S\t{n}\t{v}\n" for n, v in seg.items())
        gfa += "".join(f"L\t{a}\t{ao}\t{b}\t{bo}\t{ov}M\n" for a, ao, b, bo, ov in links)
        reads = [(circle + circle)[i * 76:i * 76 + 12000] for i in range(500)]  # >=59 per entry
        # evidence names the pieces (as an OVASM graph with a non-maximal chain would)
        evidence = {"links": [{"from": "a1q+", "to": "X1-", "reads": 7}],
                    "traversals": [{"path": "a1q+ X1- a2-", "reads": 1}],
                    "repeats": [{"repeat": "X1-", "entries": ["a1q+", "a2+"], "exits": ["X2+", "X2+"],
                                 "pairing_reads": [1, 1]},
                                {"repeat": "X2-", "entries": ["b1-", "b2-"], "exits": ["X1+", "X1+"],
                                 "pairing_reads": [1, 1]}]}
        with TemporaryDirectory() as d:
            d = Path(d)
            (d / "g.gfa").write_text(gfa)
            (d / "r.fa").write_text("".join(f">r{i}\n{r}\n" for i, r in enumerate(reads)))
            (d / "ev.json").write_text(json.dumps(evidence))
            subprocess.run([sys.executable, "-I", str(SCRIPT), "--gfa", str(d / "g.gfa"), "--reads",
                            str(d / "r.fa"), "--out", str(d / "out"), "--minimap2", MINIMAP2, "--threads", "2",
                            "--evidence", str(d / "ev.json")],
                           check=True, capture_output=True)
            links = {r.split("\t")[0]: r.split("\t") for r in (d / "out/links.tsv").read_text().splitlines()[1:]}
            self.assertIn("7", [r[5] for r in links.values()])  # a1q+ -> X1- found under its merged name
            # the repeats name a piece inside the merged repeat (X2 on one, X1 on the other): skipped
            self.assertEqual((d / "out/repeats.tsv").read_text(), "")
            summary = json.loads((d / "out/summary.json").read_text())
            self.assertEqual(summary["segments"], 5)
            self.assertEqual(summary["branching_segments"], 1)
            rows = [r.split("\t") for r in (d / "out/traversals.tsv").read_text().splitlines()[1:]]
            find = lambda a, b: [r for r in rows if a in r[0] and b in r[0]][0]
            # the merged repeat is named after its parts in walk order, whichever strand that is
            self.assertIn(find("a1p", "b1")[1], ("X1-|X2+", "X2-|X1+"))
            self.assertEqual(find("a1p", "b1")[2], "3000")
            self.assertEqual(find("a1p", "b1")[7], "supported")
            self.assertEqual(find("a2", "b2")[7], "supported")
            self.assertEqual(find("a1p", "b2")[7], "contradicted")
            self.assertEqual(find("a2", "b1")[7], "contradicted")

    def test_traversals_links_and_missed_links(self):
        rng = random.Random(7)
        s, circle = self.molecules(rng)
        gfa = "".join(f"S\t{n}\t{v}\n" for n, v in s.items())
        gfa += "".join(f"L\t{a}\t+\t{b}\t+\t0M\tev:i:5\n" for a, b in
                       (("a1", "X"), ("a2", "X"), ("X", "b1"), ("X", "b2")))
        reads = []
        for i in range(400):
            start = rng.randrange(len(circle))
            r = (circle + circle)[start:start + 12000]
            r = "".join(c if rng.random() > 0.001 else rng.choice("ACGT") for c in r)
            reads.append(r if i % 2 else r.translate(str.maketrans("ACGT", "TGCA"))[::-1])
        evidence = {"links": [], "traversals": [{"path": "a1+ X+ b1+", "reads": 9},
                                                {"path": "a1+ X+ b2+", "reads": 3}],
                    "repeats": [{"repeat": "X+", "entries": ["a1+", "a2+"], "exits": ["b1+", "b2+"],
                                 "pairing_reads": [9, 3]}]}
        with TemporaryDirectory() as d:
            d = Path(d)
            (d / "g.gfa").write_text(gfa)
            (d / "r.fa").write_text("".join(f">r{i}\n{r}\n" for i, r in enumerate(reads)))
            (d / "ev.json").write_text(json.dumps(evidence))
            subprocess.run([sys.executable, "-I", str(SCRIPT), "--gfa", str(d / "g.gfa"), "--reads",
                            str(d / "r.fa"), "--out", str(d / "out"), "--minimap2", MINIMAP2,
                            "--evidence", str(d / "ev.json"), "--threads", "2"],
                           check=True, capture_output=True)
            summary = json.loads((d / "out/summary.json").read_text())
            trav = {r.split("\t")[0]: r.split("\t") for r in
                    (d / "out/traversals.tsv").read_text().splitlines()[1:]}
            self.assertEqual(trav["a1+ X+ b1+"][7], "supported")
            self.assertEqual(trav["a2+ X+ b2+"][7], "supported")
            self.assertEqual(trav["a1+ X+ b2+"][7], "contradicted")
            self.assertEqual(trav["a2+ X+ b1+"][7], "contradicted")
            self.assertEqual(summary["link_status"], {"supported": 4})
            missed = (d / "out/missed_links.tsv").read_text()
            self.assertIn("a2->b1-", missed)  # canonical form of b1+ -> a2+
            self.assertIn("a1->b2-", missed)  # canonical form of b2+ -> a1+
            self.assertEqual(summary["missed_links"], 2)
            rep = (d / "out/repeats.tsv").read_text().splitlines()[1].split("\t")
            self.assertEqual(rep[4], "interconvertible")  # OVASM's (synthetic) call
            self.assertEqual(rep[6], "apart")             # what the independent reads show
            self.assertEqual(summary["ovasm_vs_heldout"]["repeat_calls"], {"interconvertible->apart": 1})


if __name__ == "__main__":
    unittest.main()
