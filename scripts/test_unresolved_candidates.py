import json, os, random, subprocess, tempfile, unittest
from pathlib import Path
SCRIPT=Path(__file__).with_name("assemble_unresolved_candidates.py")
BINARY=os.environ.get("OVASM_TEST_BINARY")
MINIMAP=os.environ.get("OVASM_TEST_MINIMAP2")
@unittest.skipUnless(BINARY and MINIMAP,"requires verified formal-source binary and minimap2")
class CandidateGraphTest(unittest.TestCase):
 def test_shared_gene_graph_retains_ambiguous_flanks(self):
  def sequence(seed,n): return "".join(random.Random(seed).choices("ACGT",k=n))
  with tempfile.TemporaryDirectory(dir=os.environ.get("OVASM_TEST_TMP")) as tmp:
   w=Path(tmp); gene=sequence(12,1500); refs=[]
   for name,seed,step in [("true",1,30),("high",2,10),("nuclear_same_depth",3,30)]:
    seq=sequence(seed,10000);refs.append((name,seq[:4000]+gene+seq[5500:],step))
   ref=w/"ref.fa"; ref.write_text("".join(">"+n+"\n"+s+"\n" for n,s,_ in refs));reads=w/"reads.fa"
   with reads.open("w") as f:
    for name,s,step in refs:
     for pos in range(1000,5500,step):
      overlap=max(0,min(pos+3000,5500)-max(pos,4000))
      if 1000<=overlap<2400:f.write(">"+name+"_"+str(pos)+"\n"+s[pos:pos+3000]+"\n")
   out=w/"out";subprocess.run(["python3","-I",str(SCRIPT),"--binary",BINARY,"--reads",str(reads),"--out",str(out)],check=True)
   manifest=json.loads((out/"manifest.json").read_text());self.assertIn("unresolved",manifest["classification"]);self.assertIn("not certified",manifest["classification"])
   seg=w/"segments.fa"; seg.write_text("".join(">"+c[1]+"\n"+c[2]+"\n" for c in (l.split("\t") for l in (out/"assembly.gfa").read_text().splitlines()) if c[0]=="S" and c[2]!="*"))
   result=subprocess.run([MINIMAP,"-x","asm5","-c","-N","50","-p","0.5",str(ref),str(seg)],check=True,capture_output=True,text=True)
   hits={n:[] for n,_,_ in refs}
   for line in result.stdout.splitlines():
    c=line.split("\t")
    if int(c[10])>=200 and int(c[9])/int(c[10])>=.95:hits[c[5]].append((int(c[7]),int(c[8])))
   for n in hits:
    self.assertTrue(any(a<3500 for a,b in hits[n]),n+" left flank missing")
    self.assertTrue(any(b>6000 for a,b in hits[n]),n+" right flank missing")
 def test_success_exit_without_graph_is_failure(self):
  with tempfile.TemporaryDirectory(dir=os.environ.get("OVASM_TEST_TMP")) as tmp:
   w=Path(tmp);reads=w/"r.fa";reads.write_text(">x\nACGT\n");fake=w/"fake";fake.write_text("#!/bin/sh\nexit 0\n");fake.chmod(0o755);out=w/"out"
   p=subprocess.run(["python3","-I",str(SCRIPT),"--binary",str(fake),"--reads",str(reads),"--out",str(out)],capture_output=True)
   self.assertNotEqual(p.returncode,0);self.assertIn("invalid graph output",json.loads((out/"status.json").read_text())["termination"])
 def test_existing_output_is_not_overwritten(self):
  with tempfile.TemporaryDirectory(dir=os.environ.get("OVASM_TEST_TMP")) as tmp:
   w=Path(tmp);reads=w/"r.fa";reads.write_text(">x\nACGT\n");out=w/"out";out.mkdir();marker=out/"marker";marker.write_text("keep")
   p=subprocess.run(["python3","-I",str(SCRIPT),"--binary",BINARY,"--reads",str(reads),"--out",str(out)],capture_output=True)
   self.assertNotEqual(p.returncode,0);self.assertEqual(marker.read_text(),"keep")
if __name__=="__main__": unittest.main()
