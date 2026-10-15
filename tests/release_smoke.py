"""Small offline release/installer checks against this formal repository."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import random
import shutil
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('--package-root', type=Path, required=True)
parser.add_argument('--out', type=Path, required=True)
parser.add_argument('--archive', type=Path)
parser.add_argument('--checksums', type=Path)
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
out = args.out.resolve()
out.mkdir(parents=True, exist_ok=True)
package = args.package_root.resolve()
windows = platform.system() == 'Windows'
exe = 'ovasm.exe' if windows else 'ovasm'
checks = []

def run(command, expected=0):
    result = subprocess.run([str(x) for x in command], capture_output=True, text=True)
    if (result.returncode == 0) != (expected == 0):
        raise AssertionError((command, result.returncode, result.stdout[-500:], result.stderr[-500:]))
    return result

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

if windows:
    installer = root / 'scripts/install.ps1'
    assert digest(installer) == digest(package / 'install.ps1')
    import winreg
    def user_path():
        with winreg.OpenKey(winreg.HKEY_CURRENT_USER, 'Environment') as key:
            try:
                return winreg.QueryValueEx(key, 'Path')
            except FileNotFoundError:
                return None
    before = user_path()
    prefix = ['powershell.exe', '-NoProfile', '-ExecutionPolicy', 'Bypass', '-File', installer]
    run(prefix + ['-PackageRoot', package, '-InstallDir', out / 'installed', '-NoPathUpdate'])
    assert user_path() == before
    checks.append('NoPathUpdate preserves the real user PATH')
    bad = out / 'corrupt-package'
    (bad / 'bin').mkdir(parents=True)
    shutil.copy2(package / 'FILES-SHA256.txt', bad / 'FILES-SHA256.txt')
    content = bytearray((package / 'bin' / exe).read_bytes())
    content[100] ^= 1
    (bad / 'bin' / exe).write_bytes(content)
    prior = digest(out / 'installed' / exe)
    run(prefix + ['-PackageRoot', bad, '-InstallDir', out / 'installed', '-NoPathUpdate'], expected=1)
    assert digest(out / 'installed' / exe) == prior
else:
    assert args.archive and args.checksums
    installer = root / 'scripts/install.sh'
    assert digest(installer) == digest(package / 'install.sh')
    run(['sh', '-n', installer])
    run(['sh', installer, '--archive', args.archive, '--checksums', args.checksums, '--bin-dir', out / 'installed'])
    bad = out / 'bad-checksums'
    bad.write_text('0' * 64 + '  ' + args.archive.name + '\n')
    prior = digest(out / 'installed' / exe)
    run(['sh', installer, '--archive', args.archive, '--checksums', bad, '--bin-dir', out / 'installed'], expected=1)
    assert digest(out / 'installed' / exe) == prior
    run(['sh', installer, '--version', 'invalid'], expected=1)
    checks.append('Invalid version rejected before download')
checks += ['Valid local installation succeeds', 'Corrupt package is rejected without replacing the installed executable']
installed = out / 'installed' / exe
assert digest(installed) == digest(package / 'bin' / exe)
assert run([installed, '--version']).stdout.strip() == 'ovasm 0.2.0'
run([installed, '--help'])
checks += ['Installed executable matches package hash', 'Version and help work without build tools']
# A tiny synthetic junction verifies parsing and the actual evidence path. This
# is a CLI fixture, not a biological assembly or new performance benchmark.
randomizer = random.Random(20261015)
sequence = ''.join(randomizer.choice('ACGT') for _ in range(2800))
(out / 'toy.gfa').write_text('S\tleft\t' + sequence[:1400] + '\nS\tright\t' + sequence[1400:] + '\nL\tleft\t+\tright\t+\t0M\n')
read = sequence[900:1900]
(out / 'toy.fastq').write_text(''.join(f'@fixture{i}\n{read}\n+\n' + 'I' * len(read) + '\n' for i in range(3)))
run([installed, 'evidence', '--gfa', out / 'toy.gfa', '--reads', out / 'toy.fastq', '--out', out / 'toy-evidence.json',
     '--k', '21', '--flank', '500', '--min-anchor', '200', '--bootstrap', '0', '--threads', '1'])
evidence = json.loads((out / 'toy-evidence.json').read_text())
assert evidence['reads_scanned'] == 3 and evidence['links_supported'] == 1
assert evidence['links'][0]['reads'] == 3
checks.append('Installed CLI counts all3 complete synthetic junction records')
receipt = {'status': 'passed', 'platform': platform.system(), 'formal_test_source': str(Path(__file__).resolve()),
           'formal_installer_source': str(installer), 'checks': checks, 'scientific_acceptance': False}
(out / 'receipt.json').write_text(json.dumps(receipt, indent=2))
print(json.dumps(receipt, indent=2))
