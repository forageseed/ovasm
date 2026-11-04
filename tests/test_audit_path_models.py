#!/usr/bin/env python3
"""Behavioral controls for the formal-repository path audit."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SOURCE = Path(__file__).resolve().parents[1] / 'scripts/audit_path_models.py'
spec = importlib.util.spec_from_file_location('audit_path_models', SOURCE)
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class AuditTests(unittest.TestCase):
    def test_both_strands_and_only_circular_wraps(self):
        p = ('A+', 'B+', 'C+')
        self.assertTrue(audit.contains(p, True, ('C+', 'A+', 'B+')))
        self.assertFalse(audit.contains(p, False, ('C+', 'A+', 'B+')))
        self.assertTrue(audit.contains(p, False, ('C-', 'B-', 'A-')))
        self.assertFalse(audit.contains(p, False, ('A+', 'C+', 'B+')))

    def fixture(self, directory, walks, candidates, local=()):
        d = Path(directory)
        g = d / 'graph.gfa'
        # Deliberately uncompressed graph: two repeat branch points share an internal run.
        edges = [('A+', 'R+'), ('C+', 'R+'), ('R+', 'M+'), ('M+', 'S+'), ('S+', 'B+'), ('S+', 'D+'), ('B+', 'C+'), ('D+', 'A+')]
        g.write_text(''.join(f'S\t{n}\tACGT\n' for n in ('A','B','C','D','M','R','S')) +
                     ''.join(f'L\t{a[:-1]}\t{a[-1]}\t{b[:-1]}\t{b[-1]}\t0M\n' for a,b in edges))
        e = d / 'evidence.json'
        e.write_text(json.dumps(dict(graph=str(g), links=[{'from':a,'to':b} for a,b in edges],
                                     observed_walks=[{'path':p,'reads':n} for p,n in walks],
                                     traversals=[{'path':p,'reads':n} for p,n in local])))
        covers = [dict(molecules=[dict(path=p,circular=c) for p,c in ms], weak_bridges=0,
                       topological_bridges=0, copy_deviation=0) for ms in candidates]
        l = d / 'linearization.json'
        l.write_text(json.dumps(dict(graph=str(g), components=[dict(best=covers[0],alternatives=covers[1:],
                                     candidates=len(covers),search_truncated=False,cover_size_limited=False)])))
        return g,e,l

    def test_joint_linkage_separates_equal_local_marginals(self):
        true = [('A+ R+ M+ S+ B+', False), ('C+ R+ M+ S+ D+', False)]
        false = [('A+ R+ M+ S+ D+', False), ('C+ R+ M+ S+ B+', False)]
        local = [('A+ R+ M+',10),('C+ R+ M+',10),('M+ S+ B+',10),('M+ S+ D+',10)]
        with tempfile.TemporaryDirectory() as d:
            p = self.fixture(d,[('A+ R+ M+ S+ B+',10),('C+ R+ M+ S+ D+',10)],[false,true],local)
            r = audit.audit(*p)
            c = r['components'][0]
            self.assertEqual(c['top_tied_candidates'],[1])
            self.assertTrue(c['original_excluded_from_top'])
            scores = c['ordered_candidates']
            self.assertEqual([x['matched_record_walk_occurrences']['local'] for x in scores],[40,40])
            self.assertEqual([x['matched_record_walk_occurrences']['joint'] for x in scores],[20,0])
            self.assertFalse(r['physical_topology_identified'])

    def test_coexisting_joint_configurations_stay_tied(self):
        a=[('A+ R+ M+ S+ B+',False),('C+ R+ M+ S+ D+',False)]
        b=[('A+ R+ M+ S+ D+',False),('C+ R+ M+ S+ B+',False)]
        with tempfile.TemporaryDirectory() as d:
            p=self.fixture(d,[(x,5) for x,_ in a+b],[a,b])
            self.assertEqual(audit.audit(*p)['components'][0]['top_tied_candidates'],[0,1])

    def test_no_circle_bonus_without_closing_observation(self):
        p='A+ R+ M+ S+ D+'
        with tempfile.TemporaryDirectory() as d:
            args=self.fixture(d,[(p,10)], [[(p,True)],[(p,False)]])
            self.assertEqual(audit.audit(*args)['components'][0]['top_tied_candidates'],[0,1])

    def test_observed_closing_chain_changes_evidence(self):
        p='A+ R+ M+ S+ D+'
        with tempfile.TemporaryDirectory() as d:
            args=self.fixture(d,[('S+ D+ A+ R+',10)], [[(p,False)],[(p,True)]])
            self.assertEqual(audit.audit(*args)['components'][0]['top_tied_candidates'],[1])

    def test_graph_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            args=self.fixture(d,[],[[('A+ R+ M+',False)]])
            e=json.loads(args[1].read_text()); e['links'].pop(); args[1].write_text(json.dumps(e))
            with self.assertRaisesRegex(ValueError,'links do not match'): audit.audit(*args)

    def test_recorded_graph_sequence_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            args=self.fixture(d,[],[[('A+ R+ M+',False)]])
            other=Path(d)/'other.gfa'; other.write_text(args[0].read_text().replace('ACGT','AAAA'))
            with self.assertRaisesRegex(ValueError,'recorded graph missing or differs'):
                audit.audit(other,*args[1:])

    def test_reverse_duplicate_observations_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            args=self.fixture(d,[('A+ R+ M+',10),('M- R- A-',10)],[[('A+ R+ M+',False)]])
            with self.assertRaisesRegex(ValueError,'double-count'): audit.audit(*args)

    def test_truncation_and_missing_evidence_are_preserved(self):
        with tempfile.TemporaryDirectory() as d:
            args=self.fixture(d,[],[[('A+ R+ M+',False)]])
            e=json.loads(args[1].read_text()); del e['observed_walks']; args[1].write_text(json.dumps(e))
            l=json.loads(args[2].read_text()); c=l['components'][0]
            c.update(search_truncated=True,cover_size_limited=True,candidates=99)
            args[2].write_text(json.dumps(l)); r=audit.audit(*args); c=r['components'][0]
            self.assertTrue(r['missing_walk_evidence'])
            self.assertTrue(c['search_truncated'] and c['cover_size_limited'] and c['candidate_list_limited'])
            self.assertFalse(c['has_matched_sequence_evidence'])

    def test_no_partial_ranking_when_budget_exceeded(self):
        with tempfile.TemporaryDirectory() as d:
            args=self.fixture(d,[('A+ R+ M+',10)],[[('A+ R+ M+',False)]])
            with self.assertRaisesRegex(ValueError,'budget exceeded'): audit.audit(*args,max_comparisons=1)


if __name__=='__main__':
    unittest.main(verbosity=2)
