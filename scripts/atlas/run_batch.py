"""Run `ovasm run` over a manifest, resumable, one sample after another per CPU slot.

Each sample is assembled into OUT/<sample>.tmp under /usr/bin/time -v, taskset, nice 19 and a
hard memory limit (prlimit --as), then renamed to OUT/<sample> with status.json: command, binary
sha256, exit code, wall seconds, peak RSS. A sample with status.json is skipped; a leftover
.tmp from an interrupted run is moved to OUT/.attic/ first (never deleted).

usage: run_batch.py --manifest M.tsv --out DIR --ovasm BIN [--slots 24-27,80-83 ...]
       [--threads 8] [--memory-gib 64] [--only SAMPLE ...] [--limit N] [--order small|manifest]
Each --slots value is one CPU set running one sample at a time.
"""
import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

lock = threading.Lock()


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def parse_time(path):
    out = {}
    if path.exists():
        for line in path.read_text().splitlines():
            line = line.strip()
            if line.startswith('Elapsed (wall clock) time'):
                v = line.split(': ', 1)[1].strip().split(':')
                s = 0.0
                for part in v:
                    s = s * 60 + float(part)
                out['wall_seconds'] = round(s, 1)
            elif line.startswith('Maximum resident set size'):
                out['peak_rss_gb'] = round(int(line.split(': ', 1)[1]) / 1024 ** 2, 2)
            elif line.startswith('Exit status'):
                out['exit_status'] = int(line.split(': ', 1)[1])
    return out


def log(msg):
    with lock:
        print(time.strftime('%F %T'), msg, flush=True)


def run_one(row, args, slot, binary_sha):
    out = Path(args.out)
    final = out / row['sample']
    tmp = out / (row['sample'] + '.tmp')
    # leftovers of an interrupted run (no status.json): keep them aside, never delete
    for left in (tmp, final):
        if left.exists():
            attic = out / '.attic'
            attic.mkdir(parents=True, exist_ok=True)
            left.rename(attic / f"{left.name}.{time.strftime('%Y%m%d%H%M%S')}")
    tmp.parent.mkdir(parents=True, exist_ok=True)
    sub = 'short' if row['read_type'] == 'short' else 'hifi'
    cmd = [args.ovasm, 'run', sub, '--organelle', 'both', '--threads', str(args.threads),
           '--sample', row['sample'], '--out', str(tmp / 'ovasm')]
    for r in row['reads'].split(';'):
        cmd += ['--reads', r]
    wrapped = ['/usr/bin/time', '-v', '-o', str(tmp / 'time.txt'),
               'taskset', '-c', slot, 'nice', '-n', '19',
               'prlimit', f'--as={args.memory_gib * 1024 ** 3}', '--'] + cmd
    tmp.mkdir()
    log(f'start {row["sample"]} on {slot}')
    t0 = time.time()
    with open(tmp / 'ovasm.log', 'w') as fh:
        code = subprocess.run(wrapped, stdout=fh, stderr=subprocess.STDOUT, cwd='/tmp').returncode
    status = {
        'sample': row['sample'],
        'species': row['species'],
        'source': row['source'],
        'read_type': row['read_type'],
        'command': cmd,
        'binary': args.ovasm,
        'binary_sha256': binary_sha,
        'cpu_slot': slot,
        'memory_limit_gib': args.memory_gib,
        'exit_code': code,
        'started': time.strftime('%FT%T', time.localtime(t0)),
        **parse_time(tmp / 'time.txt'),
    }
    (tmp / 'status.json').write_text(json.dumps(status, indent=1))
    tmp.rename(final)
    log(f'done {row["sample"]} exit {code} {status.get("wall_seconds")} s {status.get("peak_rss_gb")} GB')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--ovasm', required=True)
    ap.add_argument('--slots', action='append', default=None)
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--memory-gib', type=int, default=64)
    ap.add_argument('--only', action='append', default=None)
    ap.add_argument('--limit', type=int, default=None)
    ap.add_argument('--order', choices=['small', 'manifest'], default='small')
    args = ap.parse_args()
    slots = args.slots or ['24-27,80-83']
    args.ovasm = os.path.abspath(args.ovasm)
    rows = list(csv.DictReader(open(args.manifest), delimiter='\t'))
    if args.only:
        rows = [r for r in rows if r['sample'] in args.only]
    if args.order == 'small':
        rows.sort(key=lambda r: int(r.get('bytes') or 0))
    todo = [r for r in rows if not (Path(args.out) / r['sample'] / 'status.json').exists()]
    if args.limit is not None:
        todo = todo[: args.limit]
    log(f'{len(rows)} in manifest, {len(todo)} to run, slots {slots}')
    binary_sha = sha256(args.ovasm)
    queue = list(todo)

    def worker(slot):
        while True:
            with lock:
                if not queue:
                    return
                row = queue.pop(0)
            try:
                run_one(row, args, slot, binary_sha)
            except Exception as e:  # keep the batch going; the sample has no status.json
                log(f'error {row["sample"]}: {e}')

    threads = [threading.Thread(target=worker, args=(s,)) for s in slots]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    log('all done')


if __name__ == '__main__':
    sys.exit(main())
