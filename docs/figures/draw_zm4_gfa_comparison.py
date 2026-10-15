"""Reproduce the reviewed ZM4 OVASM/OATK evidence comparison from source data.

This is figure code, not an assembler, new scientific scoring or read-origin inference.
Requires matplotlib and numpy. The paired JSON retains graph hashes, display strands,
exact homology intervals and native evidence provenance; no raw reads are needed.
"""
import argparse,json,math,sys,hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import PathPatch,Patch,FancyArrowPatch
from matplotlib.path import Path as MPath
from matplotlib.transforms import Bbox

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--data',type=Path,default=Path(__file__).with_name('zm4-gfa-comparison-data.json'))
parser.add_argument('--out',type=Path,default=Path(__file__).resolve().parent)
parser.add_argument('--alignment-helper',type=Path,help='Optional external audit_panel_alignment.py, for an additional publication QA gate')
args=parser.parse_args();O=args.out.resolve();O.mkdir(parents=True,exist_ok=True)
data=json.loads(args.data.read_text(encoding='utf-8'));assert data['schema_version']==1
COL=data['colours'];INK=data['ink'];G={k:(v['nodes'],v['links']) for k,v in data['graphs'].items()};D={}
for k,nodes in data['layout'].items():
 for n,v in nodes.items():D[k,n]={**v,'points':np.asarray(v['points'],dtype=float)}
rows=data['edge_support']
assert {k:(len(G[k][0]),len(G[k][1])) for k in G}=={'ovasm':(9,12),'oatk':(12,16)}
assert all(set(G[k][0])==set(data['layout'][k]) for k in G)
assert all(len(D[k,n]['roles'])==len(D[k,n]['points'])-1 for k in G for n in G[k][0])
assert all(r['HiFi']>0 and r['ONT']>0 for r in rows)
assert data['integrity']['homology']['exact_sequence_and_orientation_verified_nodes']==20
matplotlib.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':7,'svg.fonttype':'none','pdf.fonttype':42,'axes.linewidth':.6})
def require_matplotlib_panel_alignment(fig,json_out,overlay_svg,panel_ids,column_groups,strict):
 # This case has exactly two equal, stacked graph panels. Measure final rendered
 # rectangles in points; fail rather than silently relaxing the 1.5 pt contract.
 fig.canvas.draw();rects=[ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted()).bounds for ax in axes]
 rects=[[v*72 for v in r] for r in rects];a,b=rects
 differences={'left_pt':abs(a[0]-b[0]),'right_pt':abs(a[0]+a[2]-b[0]-b[2]),'width_pt':abs(a[2]-b[2]),'height_pt':abs(a[3]-b[3])}
 assert max(differences.values())<=1.5,differences
 report={'verdict':'PASS','scope':'Two stacked comparable graph panels; final rendered rectangles','tolerance_pt':1.5,'rectangles_pt':dict(zip(panel_ids,rects)),'differences':differences}
 Path(json_out).write_text(json.dumps(report,indent=2),encoding='utf-8')
 if args.alignment_helper:
  import importlib.util
  spec=importlib.util.spec_from_file_location('optional_panel_audit',args.alignment_helper.resolve());helper=importlib.util.module_from_spec(spec);sys.modules[spec.name]=helper;spec.loader.exec_module(helper)
  helper.require_matplotlib_panel_alignment(fig,json_out=str(O/'publication-panel-alignment.json'),overlay_svg=overlay_svg,panel_ids=panel_ids,column_groups=column_groups,strict=strict)
fig=plt.figure(figsize=(7.2047244094,9.7244094488));axes=[];all_labels=[];edge_paths={}
fig.text(.035,.976,'Same backbone, two supported alternatives, one nuclear-context branch',fontsize=9,weight='bold',color=INK)
leg=[Patch(facecolor=COL[r],label=l) for r,l in [('A','A backbone'),('B','B backbone'),('R','Shared R'),('J1','J1'),('J2','J2'),('N','Nuclear context')]]
fig.legend(handles=leg,loc='upper left',bbox_to_anchor=(.025,.964),ncol=6,frameon=False,fontsize=6.5,columnspacing=1.4,handlelength=1.4)
# Node labels and edge tags use leaders; every actual adjacency is retained once.
node_labels={
 'ovasm':{'u0':(-1.82,.18),'u2':(-.23,.03),'u1__p0':(.38,.83),'u1__p1':(1.45,.70),'u1__p2':(2.04,.05),'u1__p3':(1.48,-.70),'u1__p4':(.36,-.86),'jr4_0':(.87,.42),'jr4_1':(.87,-.48)},
 'oatk':{'u619':(-1.95,-.20),'u12149':(-1.90,.53),'u71136':(-1.13,.83),'u12147':(-.40,.84),'u12146':(-.26,-.04),'u12145':(.31,.84),'u12144':(2.03,.26),'u71135':(1.53,-.65),'u12148':(.41,-.86),'u71133':(.87,.42),'u71134':(.87,-.48),'u12151':(-.93,.34)}}
alias={'u0':'A / u0','u2':'R / u2','jr4_0':'J1 / jr4_0','jr4_1':'J2 / jr4_1',**{'u1__p'+str(i):'B'+str(i+1)+' / u1__p'+str(i) for i in range(5)}}
# Link labels are placed automatically against rendered ribbons/labels below.
for panel,(k,y) in enumerate([('ovasm',.683),('oatk',.424)]):
 ax=fig.add_axes([.035,y,.93,.211]);axes.append(ax);ax.set_xlim(-2.40,2.40);ax.set_ylim(-1.05,1.05);ax.axis('off')
 fig.text(.035,y+.224,chr(97+panel),fontsize=9,weight='bold');fig.text(.072,y+.224,('OVASM repaired: 9 nodes / 12 links' if k=='ovasm' else 'OATK: 12 nodes / 16 links'),fontsize=8,weight='bold')
 ax.text(-1.15,-.28,'A\n152,685 bp',ha='center',va='center',fontsize=7,color=COL['A'])
 ax.text(2.08,-.30,'B\n137,805 bp',ha='center',va='center',fontsize=7,color=COL['B'])
 if k=='oatk':ax.text(-.93,-.035,'Nuclear context\n22/32 ONT bilateral',ha='center',va='center',fontsize=6.2,color=COL['N'])
 # Edges first; for overlap edges the short hook explicitly separates the two node ports.
 for i,(a,sa,b,sb,overlap) in enumerate(G[k][1]):
  da,db=D[k,a],D[k,b];fa=(sa==da['strand']);fb=(sb==db['strand']);pa=da['points'];pb=db['points'];p=pa[-1] if fa else pa[0];q=pb[0] if fb else pb[-1]
  ta=(pa[-1]-pa[-3 if len(pa)>2 else -2]) if fa else (pa[0]-pa[2 if len(pa)>2 else 1]);tb=(pb[min(2,len(pb)-1)]-pb[0]) if fb else (pb[-3 if len(pb)>2 else -2]-pb[-1]);ta/=max(np.linalg.norm(ta),1e-9);tb/=max(np.linalg.norm(tb),1e-9)
  step=max(.12,min(.45,np.linalg.norm(q-p)*.35));c1=p+ta*step;c2=q-tb*step
  if np.linalg.norm(q-p)<.15:
   mid=(p+q)/2;rad=mid.copy();rad[0]=rad[0] if abs(rad[0])>.15 else -.5;rad/=max(np.linalg.norm(rad),1e-9);c1=p+rad*.16;c2=q+rad*.16
  path=MPath([p,c1,c2,q],[MPath.MOVETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4]);isn='u12151' in (a,b);color=COL['N'] if isn else '#6F7982'
  ax.add_patch(PathPatch(path,fill=False,lw=.8,edgecolor=color,linestyle='--' if isn else '-',zorder=2))
  t=.5;mid=(1-t)**3*p+3*(1-t)**2*t*c1+3*(1-t)*t*t*c2+t**3*q
  eid=('V' if k=='ovasm' else 'T')+str(i+1);edge_paths[k,eid]=path
  label=ax.text(*mid,eid,fontsize=6.1,ha='center',va='center',color=color,zorder=7);all_labels.append((ax,label,mid,'edge',eid))
 for n,s in G[k][0].items():
  d=D[k,n];p=d['points'];segments=np.stack([p[:-1],p[1:]],axis=1);colors=[COL[r] for r in d['roles']];ax.add_collection(LineCollection(segments,colors=colors,linewidths=4.4,zorder=3,capstyle='round'))
  ax.plot(*p[[0,-1]].T,'o',ms=2.6,mfc='white',mec=INK,mew=.45,zorder=5)
  z=len(p)//2
  if len(p)>20:
   ax.annotate('',xy=p[min(len(p)-1,z+6)],xytext=p[max(0,z-6)],arrowprops={'arrowstyle':'-|>','color':'white','lw':.8,'mutation_scale':7,'shrinkA':0,'shrinkB':0},zorder=5)
  else:
   ax.annotate('',xy=p[1]*.85+p[0]*.15,xytext=p[0]*.85+p[1]*.15,arrowprops={'arrowstyle':'-|>','color':'white','lw':.8,'mutation_scale':7,'shrinkA':0,'shrinkB':0},zorder=5)
  name=(alias.get(n,n) if k=='ovasm' else ('J1 / '+n if n=='u71133' else 'J2 / '+n if n=='u71134' else n));sign=d['strand'].replace('-','−');name+=' '+sign
  target=p[len(p)//2] if len(p)>2 else np.mean(p,axis=0);xy=node_labels[k][n];label=ax.text(*xy,name,fontsize=6.5,ha='center',va='center',color=INK,zorder=8);all_labels.append((ax,label,target,'node',n))
 fig.text(.035,y-.006,'Arrows follow the same homology orientation; +/− labels retain the original GFA strand.',fontsize=6.2,color='#50606D')
# Resolve all text collisions with actual path geometry in physical point units, then add hairline leaders.
fig.canvas.draw();moves=[]
for ax in axes:
 texts=[x for x in all_labels if x[0] is ax];renderer=fig.canvas.get_renderer();fixed=[t.get_window_extent(renderer) for t in ax.texts if all(t is not x[1] for x in texts)]
 placed=[];bounds=ax.get_window_extent(renderer);obstacles=[]
 for coll in ax.collections:
  for segment in coll.get_segments():obstacles.append(MPath(ax.transData.transform(segment)))
 for patch in ax.patches:
  if not patch.get_fill():obstacles.append(patch.get_path().transformed(patch.get_transform()))
 # Native node names take precedence over numerical edge tags.
 for _,text,target,kind,n in sorted(texts,key=lambda z:z[3]!='node'):
  old=text.get_window_extent(renderer);origin=ax.transData.transform(text.get_position());chosen=None
  for r in [0,3,6,9,12,16,21,27,34,43,53,65]:
   for angle in [90,-90,0,180,45,135,-45,-135,*range(0,360,15)]:
    d=r*fig.dpi/72*np.array([math.cos(math.radians(angle)),math.sin(math.radians(angle))]);box=Bbox.from_extents(old.x0+d[0],old.y0+d[1],old.x1+d[0],old.y1+d[1]);pad=3.8*fig.dpi/72
    expanded=Bbox.from_extents(box.x0-pad,box.y0-pad,box.x1+pad,box.y1+pad)
    if box.x0<bounds.x0 or box.x1>bounds.x1 or box.y0<bounds.y0 or box.y1>bounds.y1:continue
    if any(expanded.overlaps(b) for b in fixed+placed):continue
    if any(path.intersects_bbox(expanded,filled=False) for path in obstacles):continue
    chosen=(d,box);break
   if chosen is not None:break
  if chosen is None:raise RuntimeError('label has no clear position '+n)
  d,box=chosen;new=ax.transData.inverted().transform(origin+d);text.set_position(new);placed.append(box)
  anchor=ax.transData.transform(target);distance=np.linalg.norm(ax.transData.transform(new)-anchor)

  moves.append({'panel':axes.index(ax),'label':text.get_text(),'kind':kind,'position':new.tolist()})
# Connect after the final label positions are known. Terminate outside each actual text
# box and route around other text boxes; a fixed centre-to-centre shrink is incorrect.
fig.canvas.draw()
for ax in axes:
 renderer=fig.canvas.get_renderer();texts=[x for x in all_labels if x[0] is ax]
 boxes={id(t):t.get_window_extent(renderer) for t in ax.texts};pad=2.2*fig.dpi/72
 for _,text,target,kind,n in texts:
  anchor=ax.transData.transform(target);box=boxes[id(text)];center=np.array([(box.x0+box.x1)/2,(box.y0+box.y1)/2])
  if kind=='role' or np.linalg.norm(center-anchor)<7*fig.dpi/72:continue
  bound=Bbox.from_extents(box.x0-pad,box.y0-pad,box.x1+pad,box.y1+pad)
  edge=np.clip(anchor,[bound.x0,bound.y0],[bound.x1,bound.y1]);other=[Bbox.from_extents(b.x0-pad,b.y0-pad,b.x1+pad,b.y1+pad) for key,b in boxes.items() if key!=id(text)]
  def clean(points):
   return not any(MPath(np.array([a,b])).intersects_bbox(c,filled=False) for a,b in zip(points,points[1:]) for c in other)
  candidates=[[anchor,edge]]
  for r in [6,12,20,30]:
   r*=fig.dpi/72
   for delta in [[r,0],[-r,0],[0,r],[0,-r]]:
    delta=np.array(delta);candidates.append([anchor,(anchor+edge)/2+delta,edge])
  route=next((q for q in candidates if clean(q)),None)
  if route is not None:
   xy=ax.transData.inverted().transform(np.array(route));ax.plot(xy[:,0],xy[:,1],lw=.35,color='#A1ABB3',zorder=1)
(O/'label-layout.json').write_text(json.dumps(moves,indent=2))
fig.text(.035,.393,'c',fontsize=9,weight='bold');fig.text(.072,.393,'Every GFA link has HiFi and ONT support; source still matters',fontsize=8,weight='bold')
# All link counts use exactly the same native evidence protocol; no assembler EC/SC depth values.
for k,x in [('ovasm',.035),('oatk',.53)]:
 fig.text(x,.375,k.upper()+' link (original GFA orientation)',fontsize=6.5,weight='bold');fig.text(x+.365,.375,'HiFi',fontsize=6.5,ha='right',weight='bold');fig.text(x+.425,.375,'ONT',fontsize=6.5,ha='right',weight='bold')
 subset=[r for r in rows if r['assembler']==k]
 for i,r in enumerate(subset):
  y=.360-i*.0107;isn='u12151' in r['from']+r['to'];color=COL['N'] if isn else INK
  # Use concise B aliases in OV table; raw identifiers retained in TSV.
  def short(s):return (alias.get(s[:-1],s[:-1]).split(' / ')[0] if k=='ovasm' else s[:-1])+s[-1].replace('-','−')
  fig.text(x,y,r['id']+'  '+short(r['from'])+' → '+short(r['to']),fontsize=6.3,color=color)
  fig.text(x+.365,y,str(r['HiFi']),ha='right',fontsize=6.3,color=color);fig.text(x+.425,y,str(r['ONT']),ha='right',fontsize=6.3,color=color)
fig.text(.035,.206,'J1 = u71133: 1,421 bp exact within jr4_0[319:1740).',fontsize=6.3,color=COL['J1'])
fig.text(.035,.194,'J2 = u71134: 1,406 bp exact within jr4_1[81:1487).',fontsize=6.3,color=COL['J2'])
fig.text(.035,.181,'These are shared local alternative connections, not missing MT sequence.',fontsize=6.3)
# Evidence decision, with strict denominator and protocol boundaries explicit.
fig.text(.035,.153,'What differs after the repair?',fontsize=8,weight='bold')
fig.text(.035,.135,'Shared: exact A/B/R backbone and both HiFi + ONT supported local alternatives (J1/J2).',fontsize=6.7)
fig.text(.035,.118,'Different: node boundaries / overlap lengths; OATK alone retains u12151 (grey dashed branch).',fontsize=6.7)
fig.text(.035,.101,'u12151: 22/32 branch-compatible ONT records have ≥500 bp target-absent nuclear context on BOTH sides.',fontsize=6.3)
fig.text(.035,.086,'Strict context audit: ≥95% whole-query coverage and identity; HiFi 44/44 nuclear-aligned, 0/44 bilateral flanks.',fontsize=6.2)
fig.text(.035,.065,'Edge counts: 29,287 recruited HiFi + 2,126 ascertained complete ONT records; records may support several links.',fontsize=6.2)
fig.text(.035,.051,'Native protocol: k21, flank 2 kb, anchor 500 bp, diagonal tolerance 30 bp, spacing tolerance 5%; no depth/EC counts.',fontsize=6.2)
fig.text(.035,.037,'Colour follows exact homology, including mixed-colour OATK nodes. R copy mapping is ambiguous; it is shown once.',fontsize=6.2)
fig.text(.035,.023,'Schematic layout; crossings are not connections. Local support does not certify physical circles or transfer direction.',fontsize=6.2)
fig.canvas.draw();require_matplotlib_panel_alignment(fig,json_out=str(O/'homology-matched.alignment.json'),overlay_svg=str(O/'homology-matched.alignment.svg'),panel_ids=['a','b'],column_groups=[['a','b']],strict=True)
fig.savefig(O/'OVASM-OATK-homology-matched.svg',dpi=600,facecolor='white')
svg_path=O/'OVASM-OATK-homology-matched.svg'
svg_path.write_text('\n'.join(line.rstrip() for line in svg_path.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
fig.savefig(O/'OVASM-OATK-homology-matched.pdf',dpi=600,facecolor='white')
fig.savefig(O/'OVASM-OATK-homology-matched.png',dpi=300,facecolor='white')
plt.close(fig)
receipt={**data['provenance'],'renderer':str(Path(__file__).resolve()),'renderer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'source_data':str(args.data.resolve()),'source_data_sha256':hashlib.sha256(args.data.read_bytes()).hexdigest(),'source_data_only':True,'raw_reads_reprocessed':False}
(O/'comparison-render-receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8');print(json.dumps({'status':'completed','renderer':receipt['renderer'],'out':str(O),'nodes':receipt['drawn_nodes'],'links':receipt['drawn_links']}))
