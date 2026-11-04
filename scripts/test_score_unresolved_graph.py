import unittest
from score_unresolved_graph import score
class ScoreTests(unittest.TestCase):
 def test_duplicate_secondary_and_reverse(self):
  def row(q,strand,a,b,match=1000,tags='tp:A:P\tcg:Z:1000M'):
   return f'{q}\t1000\t0\t1000\t{strand}\ts\t3000\t{a}\t{b}\t{match}\t1000\t60\t{tags}'
  lines=[row('x','+',0,1000),row('x','-',0,1000),row('y','-',1000,2000),row('z','+',2000,3000,900),row('w','+',2000,3000,tags='tp:A:S\tcg:Z:1000M')]
  s=score(lines,{'s':3000})['s'];self.assertEqual(s['unique_reads'],2);self.assertEqual(s['covered_span_bp'],2000);self.assertEqual(s['minimum_span_depth'],0);self.assertEqual(s['reads_touching_left_window'],1);self.assertEqual(s['reads_touching_right_window'],0)
if __name__=='__main__':unittest.main()
