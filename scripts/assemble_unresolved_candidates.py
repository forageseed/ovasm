#!/usr/bin/env python3
"""Build an evidence-only graph from unclassified candidates; does not establish source."""
import argparse,hashlib,json,os,resource,shutil,signal,subprocess,time
from pathlib import Path

def digest(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for data in iter(lambda:f.read(1048576),b""): h.update(data)
 return h.hexdigest()

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument("--binary",required=True,type=Path);p.add_argument("--reads",required=True,type=Path);p.add_argument("--out",required=True,type=Path)
 p.add_argument("--threads",type=int,default=4);p.add_argument("--k",type=int,default=101);p.add_argument("--min-count",type=int,default=3)
 p.add_argument("--max-input-bytes",type=int,default=1073741824);p.add_argument("--memory-gib",type=int,default=16);p.add_argument("--timeout-seconds",type=int,default=1800)
 a=p.parse_args()
 if not (1<=a.threads<=8 and 3<=a.k<=1001 and a.k%2==1 and a.min_count>=3 and 1<=a.memory_gib<=64 and 1<=a.timeout_seconds<=7200 and 1<=a.max_input_bytes<=1073741824): p.error("resource bounds or graph parameters out of range")
 if hasattr(os,"sched_getaffinity"):
  cpus=sorted(os.sched_getaffinity(0))[:8];os.sched_setaffinity(0,cpus)
 else: cpus=None
 if not a.reads.is_file() or a.reads.stat().st_size>a.max_input_bytes: raise ValueError("missing input or input byte limit exceeded")
 if not a.binary.is_file():raise ValueError("binary missing")
 if a.out.exists():raise FileExistsError(a.out)
 if shutil.disk_usage(a.out.parent).free<30*1024**3:raise RuntimeError("disk reserve below 30 GiB")
 a.out.mkdir(); cmd=[str(a.binary),"assemble","graph","--reads",str(a.reads),"--k",str(a.k),"--min-count",str(a.min_count),"--threads",str(a.threads),"--out",str(a.out/"assembly.gfa"),"--out-json",str(a.out/"assembly.json")]
 manifest={"classification":"unresolved evidence only; not certified mitochondrial, not inserted into recruitment bait","script_sha256":digest(Path(__file__)),"input":str(a.reads.resolve()),"input_sha256":digest(a.reads),"input_bytes":a.reads.stat().st_size,"binary":str(a.binary.resolve()),"binary_sha256":digest(a.binary),"command":cmd,"timeout_seconds":a.timeout_seconds,"memory_gib":a.memory_gib,"disk_reserve_gib":30,"cpu_affinity":cpus,"note":"short read length and shared repeats may prevent unique graph reconstruction; no nuclear-origin validation"}
 (a.out/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
 def limits():
  os.nice(19);resource.setrlimit(resource.RLIMIT_AS,(a.memory_gib*1024**3,)*2);os.setsid()
 started=time.monotonic(); reason=None
 with (a.out/"run.log").open("w") as log:
  proc=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,preexec_fn=limits)
  while proc.poll() is None:
   if time.monotonic()-started>a.timeout_seconds:reason="timeout"
   elif shutil.disk_usage(a.out).free<30*1024**3:reason="disk reserve"
   if reason:
    os.killpg(proc.pid,signal.SIGTERM)
    try:proc.wait(timeout=10)
    except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL)
    break
   time.sleep(1)
  rc=proc.wait()
 if not reason and rc==0:
  graph=a.out/"assembly.gfa"; report=a.out/"assembly.json"
  try:
   if not graph.is_file() or not graph.stat().st_size: raise ValueError("missing or empty GFA")
   if not any(len(c)>=3 and c[0]=="S" and c[2] not in ("", "*") for c in (line.rstrip("\n").split("\t") for line in graph.open())): raise ValueError("GFA has no sequence segments")
   if not report.is_file():raise ValueError("missing assembly JSON")
   json.loads(report.read_text())
  except (OSError,ValueError) as exc:reason="invalid graph output: "+str(exc)
 changed=digest(a.reads)!=manifest["input_sha256"] or digest(a.binary)!=manifest["binary_sha256"]
 if changed:reason="input or executable changed during execution"
 result={"exit_code":rc,"termination":reason,"wall_seconds":round(time.monotonic()-started,2),"graph_sha256":digest(a.out/"assembly.gfa") if (a.out/"assembly.gfa").is_file() else None}
 (a.out/"status.json").write_text(json.dumps(result,indent=2)+"\n")
 if rc or reason:raise SystemExit(1)
if __name__=="__main__":main()
