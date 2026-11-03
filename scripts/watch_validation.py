"""Wait for known validation processes, then compare all available results."""
import subprocess
import sys
import time
from pathlib import Path
root, refs, mm, output, *pids = sys.argv[1:]
# Bound lifetime even if an unrelated process eventually reuses a PID.
for attempt in range(1440):
    live = []
    for pid in pids:
        try:
            cmd = Path(f"/proc/{pid}/cmdline").read_bytes()
        except FileNotFoundError:
            continue
        if b"ovasm-d0c4c23" in cmd and root.encode() in cmd:
            live.append(pid)
    if not live: break
    time.sleep(60)
else:
    raise TimeoutError("validation still running after 24 hours")
subprocess.run([sys.executable,"-I",str(Path(__file__).with_name("compare_validation.py")), root,refs,mm,output],check=True)
