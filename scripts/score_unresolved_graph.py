#!/usr/bin/env python3
import argparse,hashlib,json,shutil,subprocess,tempfile
from pathlib import Path
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def score(lines,lengths,exclude=frozenset()):
 best={}
 for line in lines:
  c=line.split('\t')
  if len(c)<12:raise ValueError('bad PAF')
  q,ql,qs,qe,strand,t,tl,ts,te,matches,alen,mq=c[:12]
  if q in exclude:continue
  if t not in lengths or int(tl)!=lengths[t] or not 0<=int(ts)<int(te)<=int(tl):raise ValueError('bad reference or coordinates')
  if 'tp:A:P' not in c[12:] or not any(v.startswith('cg:Z:') for v in c[12:]) or int(alen)<1000 or int(matches)/int(alen)<.95:continue
  key=(q,t);v=(int(matches),int(alen),-int(ts),int(ts),int(te))
  if key not in best or v[:3]>best[key][:3]:best[key]=v
 out={}
 for t,n in lengths.items():
  e={0:0,n:0};left=right=total=0
  for (q,target),v in best.items():
   if target!=t:continue
   a,b=v[3:];e[a]=e.get(a,0)+1;e[b]=e.get(b,0)-1;total+=1;left+=a<min(500,n);right+=b>max(0,n-500)
  prev=dep=area=covered=0;minimum=None
  for x in sorted(e):
   z=x-prev
   if z:area+=z*dep;covered+=z*(dep>0);minimum=dep if minimum is None else min(minimum,dep)
   dep+=e[x];prev=x
  out[t]={'length':n,'unique_reads':total,'covered_span_bp':covered,'mean_span_depth':area/n,'minimum_span_depth':minimum,'reads_touching_left_window':left,'reads_touching_right_window':right}
 return out
def main():
 p=argparse.ArgumentParser()
 for n in ('gfa','reads','minimap2','out'):p.add_argument('--'+n,required=True,type=Path)
 p.add_argument('--exclude-reads',type=Path,help='construction FASTQ; report additional support excluding these IDs');p.add_argument('--threads',type=int,default=4);p.add_argument('--timeout',type=int,default=600);a=p.parse_args()
 if not 1<=a.threads<=8 or not 1<=a.timeout<=3600:p.error('invalid limits')
 files=(a.gfa,a.reads,a.minimap2)+((a.exclude_reads,) if a.exclude_reads else ())
 for f in files:
  if not f.is_file():raise FileNotFoundError(f)
 if a.gfa.stat().st_size>1024**3 or a.reads.stat().st_size>1024**3:raise ValueError('input exceeds 1 GiB')
 if a.out.exists():raise FileExistsError(a.out)
 if shutil.disk_usage(a.out.parent).free<30*1024**3:raise RuntimeError('disk reserve below 30 GiB')
 seq={}
 for line in a.gfa.read_text().splitlines():
  c=line.split('\t')
  if c[0]=='S':
   if len(c)<3 or c[2]=='*' or not c[2] or c[1] in seq:raise ValueError('bad segment')
   seq[c[1]]=c[2]
 if not seq:raise ValueError('empty graph')
 before={str(f):sha(f) for f in files}
 exclude=set()
 if a.exclude_reads:
  with a.exclude_reads.open() as f:
   while True:
    h=f.readline()
    if not h:break
    record_seq=f.readline();plus=f.readline();qual=f.readline()
    if not h.startswith('@') or not plus.startswith('+') or not record_seq or not qual:raise ValueError('invalid four-line FASTQ')
    exclude.add(h[1:].split()[0])
 with tempfile.TemporaryDirectory(dir=a.out.parent) as tmp:
  fa=Path(tmp)/'segments.fa';fa.write_text(''.join('>'+k+'\n'+v+'\n' for k,v in seq.items()))
  cmd=[str(a.minimap2),'-x','map-hifi','-c','-t',str(a.threads),str(fa),str(a.reads)]
  x=subprocess.run(cmd,check=True,capture_output=True,text=True,timeout=a.timeout)
  lengths={k:len(v) for k,v in seq.items()}
  segments=score(x.stdout.splitlines(),lengths)
  nonconstruction=score(x.stdout.splitlines(),lengths,exclude) if a.exclude_reads else None
 if before!={str(f):sha(f) for f in files}:raise RuntimeError('input changed')
 report={'classification':'source unresolved; construction reads are not independent validation','semantics':'PAF span depth includes deletions; endpoint touches do not prove graph-edge or circular junction support','coordinates':'zero-based half-open','filters':{'primary_only':True,'min_alignment_bp':1000,'min_identity':.95,'cigar_required':True,'endpoint_window_bp':500},'limitations':'secondary mappings excluded; primary mappings may be ambiguous in repeats; retain one best hit per read per segment','inputs_sha256':before,'command':cmd,'segments':segments,'excluding_construction_reads':nonconstruction,'exclude_read_ids':len(exclude),'independence_note':'different reads may be from same library; not an independent sample or origin validation'}
 with a.out.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
if __name__=='__main__':main()
