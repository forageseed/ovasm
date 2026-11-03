from pathlib import Path
import subprocess,json,re,collections,time
import argparse
parser=argparse.ArgumentParser(description="Compare completed pan reconstructions against diagnostic FASTA references; not physical-topology truth")
parser.add_argument('--run',action='append',required=True,help='LABEL=RECONSTRUCTION_DIRECTORY')
parser.add_argument('--reference',action='append',required=True,help='ORGANELLE=FASTA')
parser.add_argument('--gaps-json',type=Path)
parser.add_argument('--out',type=Path,required=True)
parser.add_argument('--minimap2',default='minimap2')
parser.add_argument('--threads',type=int,default=4)
args=parser.parse_args()
def mapping(specs):
 result={}
 for spec in specs:
  key,path=spec.split('=',1)
  if not key or not all(c.isalnum() or c in '_.-' for c in key) or key in result:parser.error('labels must be unique and safe')
  result[key]=Path(path)
 return result
runs=mapping(args.run);references=mapping(args.reference);out=args.out
if out.exists():parser.error('output directory already exists')
if args.threads<1:parser.error('threads must be positive')
out.mkdir(parents=True)
def fasta(path):
 d={};name=None
 for line in path.read_text().splitlines():
  if line.startswith('>'):name=line[1:].split()[0];d[name]=''
  elif name:d[name]+=line.strip()
 return d
def union(v):
 end=-1;total=0
 for a,b in sorted(v):total+=max(0,b-max(a,end));end=max(end,b)
 return total
results={}
for run,d in runs.items():
 if not (d/'reconstruction.json').exists():results[run]={'pending':True};continue
 rep=json.loads((d/'reconstruction.json').read_text());results[run]={}
 for org in ['mitochondrion','plastid']:
  ref=references.get(org);query=d/(org+'.fasta')
  if ref is None:continue
  if not query.exists():results[run][org]={'missing':True};continue
  paf=out/f'{run}.{org}.paf'
  with paf.open('w') as f,(out/f'{run}.{org}.align.log').open('w') as e:
   subprocess.run([args.minimap2,'-x','asm5','--secondary=yes','-N','20','-c','--eqx','-t',str(args.threads),str(ref),str(query)],stdout=f,stderr=e,check=True)
  covered=collections.defaultdict(list);exact=collections.defaultdict(list);matches=0;columns=0
  for line in paf.read_text().splitlines():
   f=line.split('\t')
   if int(f[10])<500:continue
   covered[f[5]].append((int(f[7]),int(f[8])));matches+=int(f[9]);columns+=int(f[10]);pos=int(f[7]);cigar=next((x[5:] for x in f[12:] if x.startswith('cg:Z:')),'')
   for n,op in re.findall(r'(\d+)([MIDNSHP=X])',cigar):
    n=int(n)
    if op=='=':exact[f[5]].append((pos,pos+n))
    if op in 'MDN=X':pos+=n
  rseq=fasta(ref);qseq=fasta(query)
  result={'reference':str(ref),'reference_records':{n:{'length':len(s),'aligned_bases':union(covered[n]),'exact_bases':union(exact[n])} for n,s in rseq.items()},'query_lengths':{n:len(s) for n,s in qseq.items()},'alignment_weighted_identity':matches/max(columns,1),'final':rep['organelles'][org].get('final')}
  if org=='mitochondrion' and args.gaps_json:
   gaps=json.loads(args.gaps_json.read_text())['existing_HiFi_candidate_sequence_comparison']
   result['historically_missing_intervals']=[{'record':n,'start':a,'end':b,'exact_bases':union([(max(a,x),min(b,y)) for x,y in exact[n] if min(b,y)>max(a,x)])} for n,entry in gaps.items() for a,b in entry['unmatched_intervals']]
  results[run][org]=result
(out/'quality-comparison.json').write_text(json.dumps({'runs':results,'caveat':'HiFi assembly is diagnostic reference, not same-individual truth. Secondary alignments included; coverage unions avoid double-counting; identity is alignment-weighted.'},indent=2))
