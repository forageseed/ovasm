"""One full-input Rust multichrom comparison against retained, previously seen evidence."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import time
from score_independent_multichrom import score


def main():
    ap=argparse.ArgumentParser(__doc__)
    ap.add_argument('--binary',type=Path,required=True)
    ap.add_argument('--prior',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--minimap2',type=Path,required=True)
    a=ap.parse_args()
    root=Path(__file__).resolve().parents[2]
    a.out=a.out.resolve(); a.binary=a.binary.resolve(); a.prior=a.prior.resolve()
    assert a.binary.is_relative_to(root/'target')
    a.out.mkdir(parents=True,exist_ok=False)
    prior=json.loads((a.prior/'state.json').read_text())
    reads=Path(prior['input']); before=reads.stat()
    old=a.prior/'published/normalized'
    old_report=json.loads((old/'assembly.json').read_text())
    assert (old_report['reads'],old_report['bases'])==(1044826,16825074678)
    assert (old_report['k'],old_report['s'],old_report['min_count'])==(1001,31,80)
    assert a.minimap2.is_file()
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    source={str(p.relative_to(root)):sha(p)for p in [root/'Cargo.toml',root/'Cargo.lock',*sorted((root/'src').glob('*.rs'))]}
    manifest={'run_id':'ovasm-rust-multichrom-Schoenoplectus-full-v021-v1',
        'formal_source':str(root),'source_hashes':source,'binary':str(a.binary),'binary_sha256':sha(a.binary),
        'input':str(reads),'input_size':before.st_size,'input_mtime_ns':before.st_mtime_ns,
        'expected_reads':1044826,'expected_bases':16825074678,
        'k':1001,'s':31,'min_count':80,'threads':4,
        'reference':str(a.prior/'truth.fasta'),'reference_scope':'post-assembly scoring only; previously seen 8-record reference, not a blind test',
        'retained_control_source_revision':prior['source_revision'],
        'retained_control_run_manifest':prior['run_manifest_id'],
        'control_graph_sha256':sha(old/'assembly.gfa'),'control_fasta_sha256':sha(old/'assembly.fasta'),
        'rss_limit_gib':96,'stage_timeout_seconds':10800}
    manifest['manifest_id']='sha256:'+hashlib.sha256(json.dumps(manifest,sort_keys=True).encode()).hexdigest()
    (a.out/'manifest.json').write_text(json.dumps(manifest,indent=2))
    state={'status':'running','controller_pid':os.getpid(),'manifest_id':manifest['manifest_id'],
           'steps':[],'physical_molecule_count':None,'physical_topology_identified':False}
    def save():
        state['heartbeat_unix']=time.time()
        tmp=a.out/'state.tmp';tmp.write_text(json.dumps(state,indent=2));tmp.replace(a.out/'state.json')
    def rss_group(pgid):
        total=0
        for path in Path('/proc').glob('[0-9]*/stat'):
            try:
                fields=path.read_text().rsplit(')',1)[1].split()
                if int(fields[2])==pgid:total+=int(fields[21])*os.sysconf('SC_PAGE_SIZE')
            except (OSError,ValueError,IndexError):pass
        return total
    def run(name,argv):
        log=a.out/(name+'.log'); timing=a.out/(name+'.time')
        started=time.monotonic()
        step={'name':name,'argv':list(map(str,argv)),'status':'running'}
        state['steps'].append(step);state['stage']=name;save()
        with log.open('w')as stream:
            process=subprocess.Popen(['/usr/bin/time','-v','-o',str(timing),*map(str,argv)],stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
            step['provider_pid']=process.pid;step['process_group']=process.pid;save()
            peak=0;reason=None
            while process.poll()is None:
                peak=max(peak,rss_group(process.pid))
                step['elapsed_seconds']=time.monotonic()-started;step['peak_group_rss_bytes_sampled']=peak
                save()
                if peak>96*1024**3:reason='96GiB sampled process-group RSS limit'
                if step['elapsed_seconds']>10800:reason='3h stage wall limit'
                if reason:
                    assert os.getpgid(process.pid)==process.pid
                    os.killpg(process.pid,signal.SIGTERM)
                    try:process.wait(timeout=10)
                    except subprocess.TimeoutExpired:os.killpg(process.pid,signal.SIGKILL);process.wait()
                    break
                time.sleep(1)
            step.update(exit_code=process.returncode,status='completed'if process.returncode==0 and not reason else 'failed',
                        wall_seconds=time.monotonic()-started,resource_stop=reason)
        text=timing.read_text()if timing.exists()else ''
        m=re.search(r'Maximum resident set size \(kbytes\):\s+(\d+)',text)
        step['gnu_max_rss_kib']=int(m[1])if m else None
        save()
        if step['status']!='completed':raise RuntimeError(f'{name}: {reason or process.returncode}; {log.read_text()[-1000:]}')
    save()
    try:
        run('rust-full-assembly',[a.binary,'assemble','--multichrom-graph','--reads',reads,'--k',1001,'--s',31,'--min-count',80,'--threads',4,'--out',a.out/'raw.gfa','--out-json',a.out/'assembly.json'])
        report=json.loads((a.out/'assembly.json').read_text())
        assert report['graph_engine']=='ovasm-rust-hpc'and not report['uses_syncasm']
        assert(report['reads'],report['bases'])==(1044826,16825074678)
        assert(reads.stat().st_size,reads.stat().st_mtime_ns)==(before.st_size,before.st_mtime_ns)
        state['assembly']={k:v for k,v in report.items()if k not in ['compressed_backbone','consensus']}
        run('components',[a.binary,'components','--gfa',a.out/'raw.gfa','--organelle','mitochondrion','--reject-other-only','--out',a.out/'mt.gfa','--out-json',a.out/'components.json'])
        run('primary',[a.binary,'primary','--gfa',a.out/'mt.gfa','--out',a.out/'all-candidates.fasta','--out-json',a.out/'primary.json'])
        run('read-evidence',[a.binary,'evidence','--gfa',a.out/'mt.gfa','--reads',reads,'--bootstrap',0,'--threads',4,'--out',a.out/'evidence.json'])
        state['stage']='post-assembly-reference-scoring';save()
        for label,candidate in [('rust',a.out/'all-candidates.fasta'),('retained-control',old/'assembly.fasta')]:
            state[label+'_score']=score(a.minimap2,a.prior/'truth.fasta',candidate,set(),a.out/(label+'-score'))
            save()
        state['status']='completed'
    except Exception as e:
        state.update(status='failed',error=repr(e))
    state['source_unchanged']=all(sha(root/p)==digest for p,digest in source.items())
    state['binary_unchanged']=sha(a.binary)==manifest['binary_sha256']
    state['stage']='terminal';save()
    print(json.dumps({k:v for k,v in state.items()if not k.endswith('_score')},indent=2))
    if state['status']!='completed':raise SystemExit(1)


if __name__=='__main__':main()
