#!/usr/bin/env python3
"""Serial pre-registered pilot. Default dry-run; --execute runs existing verified subsets."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import signal
import subprocess
import time


def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(1048576),b""):h.update(block)
    return h.hexdigest()


def group_rss(pgid):
    total=0
    for p in Path("/proc").glob("[0-9]*"):
        try:
            fields=(p/"stat").read_text().rsplit(")",1)[1].split()
            if int(fields[2]) == pgid:
                total += int((p/"statm").read_text().split()[1])*os.sysconf("SC_PAGE_SIZE")
        except (OSError,ValueError,IndexError):pass
    return total


def stop_group(pid):
    try:os.killpg(pid,signal.SIGKILL)
    except ProcessLookupError:pass


def execute(command, dest, cfg, is_ovasm):
    def setup():
        os.setsid();os.nice(19);os.sched_setaffinity(0,cfg["cpus"])
        if is_ovasm:resource.setrlimit(resource.RLIMIT_AS,(cfg["memory_bytes"],cfg["memory_bytes"]))
    env=os.environ.copy();env["PATH"]=str(Path(cfg["oatk"]).parent)+os.pathsep+env.get("PATH","")
    start=time.monotonic();peak=0;reason=None
    with (dest/"stdout.log").open("wb") as out,(dest/"stderr.log").open("wb") as err:
        proc=subprocess.Popen(["/usr/bin/time","-v","-o",str(dest/"time.txt"),*command],cwd=dest,stdout=out,stderr=err,env=env,preexec_fn=setup)
        try:
            while proc.poll() is None:
                peak=max(peak,group_rss(proc.pid))
                if peak>cfg["memory_bytes"]:reason="rss_limit"
                elif time.monotonic()-start>cfg["timeout_seconds"]:reason="timeout"
                elif shutil.disk_usage(dest).free<cfg["reserve_bytes"]:reason="disk_reserve"
                if reason:stop_group(proc.pid);break
                time.sleep(1)
            code=proc.wait()
        finally:stop_group(proc.pid)
    return dict(exit_code=code,termination=reason,wall_seconds=time.monotonic()-start,
                sampled_process_group_peak_rss_bytes=peak,memory_note="1-second RSS sum, shared pages may be counted twice; not cgroup hard limit")


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("manifest",type=Path);ap.add_argument("subsets",type=Path)
    ap.add_argument("--out",type=Path,required=True);ap.add_argument("--execute",action="store_true")
    args=ap.parse_args();cfg=json.loads(args.manifest.read_text());subpath=args.subsets.resolve()
    data=json.loads(subpath.read_text());root=args.out.resolve()
    if data.get("status")!="complete" or (subpath.parent/"INCOMPLETE").exists():raise ValueError("incomplete subsampling")
    if root.exists():raise FileExistsError(root)
    expected=cfg["reads"]
    if [str(Path(x["path"]).resolve()) for x in data["sources"]] != expected:raise ValueError("input list mismatch")
    if data["seeds"] != cfg["seeds"] or set(data["fractions"]) != set(cfg["fractions"]):raise ValueError("subset matrix mismatch")
    for key in ["ovasm","oatk","nhmmscan","mt_hmm","pt_hmm"]:
        if not Path(cfg[key]).is_file():raise FileNotFoundError(cfg[key])
    # Freeze expected binary hashes in the reviewed manifest before allowing execution.
    hashes={k:digest(cfg[k]) for k in ["ovasm","oatk","nhmmscan","mt_hmm","pt_hmm"]}
    if args.execute and hashes != cfg.get("sha256"):raise ValueError("binary/database hash mismatch: freeze reviewed hashes first")
    for key in ["mt_hmm", "pt_hmm"]:
        for suffix in [".h3f", ".h3i", ".h3m", ".h3p"]:
            db = Path(cfg[key] + suffix)
            if not db.is_file(): raise FileNotFoundError(db)
    versions = {}
    for key in ["ovasm", "oatk"]:
        version = subprocess.run([cfg[key], "--version"], capture_output=True, text=True, timeout=15, check=True)
        versions[key] = version.stdout + version.stderr
    plan=[]
    for subset in data["subsets"]:
        fq=(subpath.parent/subset["path"]).resolve()
        if not fq.is_file() or fq.stat().st_size!=subset["compressed_bytes"]:raise ValueError("missing/changed subset")
        if not subset["bases"]:raise ValueError("empty subset")
        for track in ["ovasm-default","oatk-c30","oatk-c3"]:
            dest=root/f"seed-{subset['seed']}.f-{subset['fraction']}.{track}"
            if track=="ovasm-default":cmd=[cfg["ovasm"],"run","hifi","--reads",str(fq),"--organelle","both","--threads","8","--out",str(dest/"assembly")]
            else:cmd=[cfg["oatk"],"-t","8","-c","30" if track=="oatk-c30" else "3","--nhmmscan",cfg["nhmmscan"],"-m",cfg["mt_hmm"],"-p",cfg["pt_hmm"],"-o",str(dest/"assembly"),str(fq)]
            plan.append(dict(track=track,subset=subset,subset_path=str(fq),directory=str(dest),command=cmd))
    record=dict(config=cfg,versions=versions,actual_sha256=hashes,subset_manifest=data,subset_manifest_sha256=digest(subpath),jobs=plan)
    if not args.execute:print(json.dumps(record,indent=2));return
    if shutil.disk_usage(root.parent).free<cfg["reserve_bytes"]:raise RuntimeError("disk reserve unavailable")
    root.mkdir();(root/"plan.json").write_text(json.dumps(record,indent=2)+"\n")
    # Verify compressed hashes immediately before each run; hashes also audit identical paired input.
    input_hashes = {job["subset_path"]: digest(job["subset_path"]) for job in plan}
    for job in plan:
        if shutil.disk_usage(root).free<cfg["reserve_bytes"]:raise RuntimeError("disk reserve unavailable; prior failures retained")
        dest=Path(job["directory"]);dest.mkdir();job["input_compressed_sha256"]=digest(job["subset_path"])
        if job["input_compressed_sha256"] != input_hashes[job["subset_path"]]:raise ValueError("subset changed between paired runs")
        job["result"]=execute(job["command"],dest,cfg,job["track"]=="ovasm-default")
        job["outputs"]=[dict(path=str(p.relative_to(dest)),bytes=p.stat().st_size) for p in dest.rglob("*") if p.is_file() and p.suffix.lower() in [".gfa",".fa",".fasta"] and p.stat().st_size]
        job["status"]="executed" if job["result"]["exit_code"]==0 else "failed"
        job["quality_note"]="Nonempty output and exit 0 are not assembly correctness; score independently."
        (dest/"status.json").write_text(json.dumps(job,indent=2)+"\n")
    (root/"results.json").write_text(json.dumps(record,indent=2)+"\n")

if __name__=="__main__":main()
