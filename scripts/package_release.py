#!/usr/bin/env python3
"""Package a locally built release executable, its installer and license notices."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import struct
import subprocess
import tarfile
import tomllib
import zipfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--binary', type=Path, required=True)
parser.add_argument('--out', type=Path, required=True)
parser.add_argument('--allow-dirty', action='store_true', help='Development install checks only; not for publishing')
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
version = tomllib.loads((root / 'Cargo.toml').read_text(encoding='utf-8'))['package']['version']
system = platform.system()
if system not in ('Windows', 'Linux') or platform.machine().lower() not in ('amd64', 'x86_64'):
    raise SystemExit('Only Linux x86_64 and Windows x86_64 packages are supported')
target = 'windows-x86_64' if system == 'Windows' else 'linux-x86_64'
exe = 'ovasm.exe' if system == 'Windows' else 'ovasm'
binary = args.binary.resolve(strict=True)
assert subprocess.check_output([str(binary), '--version'], text=True).strip() == 'ovasm ' + version
commit = subprocess.check_output(['git', '-C', str(root), 'rev-parse', 'HEAD'], text=True).strip()
dirty = subprocess.run(['git', '-C', str(root), 'diff', '--quiet']).returncode != 0
dirty |= subprocess.run(['git', '-C', str(root), 'diff', '--cached', '--quiet']).returncode != 0
if dirty and not args.allow_dirty:
    raise SystemExit('Commit reviewed source before producing publication assets')
out = args.out.resolve()
out.mkdir(parents=True, exist_ok=True)
name = f'ovasm-v{version}-{target}'
stage = out / name
if stage.exists():
    raise SystemExit(f'Output package directory already exists: {stage}')
(stage / 'bin').mkdir(parents=True)
shutil.copy2(binary, stage / 'bin' / exe)
if system == 'Linux':
    (stage / 'bin' / exe).chmod(0o755)
shutil.copy2(root / 'LICENSE', stage / 'LICENSE')
shutil.copy2(root / 'docs/installation.md', stage / 'INSTALL.md')
for script in (('install.cmd', 'install.ps1') if system == 'Windows' else ('install.sh',)):
    shutil.copy2(root / 'scripts' / script, stage / script)
info = {'version': version, 'target': target, 'source_commit': commit,
        'source_directory': str(root), 'tracked_source_dirty': dirty,
        'binary_sha256': hashlib.sha256(binary.read_bytes()).hexdigest(),
        'native_multichrom_kernel': system == 'Linux',
        'windows_static_crt': system == 'Windows'}
if system == 'Windows':
    image = binary.read_bytes()
    pe = struct.unpack_from('<I', image, 60)[0]
    assert image[pe:pe+4] == b'PE\0\0' and struct.unpack_from('<H', image, pe+4)[0] == 0x8664
    count = struct.unpack_from('<H', image, pe+6)[0]
    optional = pe + 24
    assert struct.unpack_from('<H', image, optional)[0] == 0x20b
    sections = optional + struct.unpack_from('<H', image, pe+20)[0]
    def file_offset(rva):
        for i in range(count):
            size, address, raw_size, raw = struct.unpack_from('<IIII', image, sections + 40*i + 8)
            if address <= rva < address + max(size, raw_size):
                return raw + rva - address
        raise ValueError('PE RVA outside sections')
    descriptor = file_offset(struct.unpack_from('<I', image, optional+120)[0])
    dependencies = []
    while any(image[descriptor:descriptor+20]):
        name = file_offset(struct.unpack_from('<I', image, descriptor+12)[0])
        dependencies.append(image[name:image.index(b'\0', name)].decode('ascii'))
        descriptor += 20
    assert not any(n.lower().startswith(('vcruntime', 'msvcp')) for n in dependencies), dependencies
    info['imported_windows_libraries'] = dependencies
else:
    symbols = subprocess.check_output(['objdump', '-T', str(binary)], text=True)
    required = max(tuple(map(int, v.split('.'))) for v in re.findall(r'GLIBC_(\d+\.\d+)', symbols))
    info['minimum_glibc'] = '.'.join(map(str, required))
    assert required <= (2, 34), 'README compatibility contract must match the built binary'
(stage / 'BUILDINFO.json').write_text(json.dumps(info, indent=2), encoding='utf-8')
# Preserve the license texts of compiled dependencies without downloading packages.
cargo = shutil.which('cargo') or str(Path.home() / '.cargo/bin/cargo')
metadata = json.loads(subprocess.check_output(
    [cargo, 'metadata', '--manifest-path', str(root / 'Cargo.toml'), '--locked', '--offline',
     '--filter-platform', 'x86_64-pc-windows-msvc' if system == 'Windows' else 'x86_64-unknown-linux-gnu',
     '--format-version', '1'],
    text=True))
notices = []
for package in metadata['packages']:
    if package['name'] == 'ovasm':
        continue
    notices.append(f"{package['name']} {package['version']}: {package.get('license') or 'see supplied license'}")
    source = Path(package['manifest_path']).parent
    dest = stage / 'licenses' / (package['name'] + '-' + package['version'])
    for p in source.iterdir():
        if p.is_file() and (p.name.upper().startswith(('LICENSE', 'LICENCE', 'COPYING', 'NOTICE'))):
            dest.mkdir(parents=True, exist_ok=True)
            shutil.copy2(p, dest / p.name)
for p in (root / 'vendor').rglob('*'):
    if p.is_file() and p.name.upper().startswith(('LICENSE', 'LICENCE', 'COPYING', 'NOTICE')):
        dest = stage / 'licenses' / 'vendor' / p.relative_to(root / 'vendor')
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, dest)
(stage / 'THIRD-PARTY-NOTICES.txt').write_text('\n'.join(notices) + '\n', encoding='utf-8')
files = sorted(p for p in stage.rglob('*') if p.is_file())
manifest = ''.join(hashlib.sha256(p.read_bytes()).hexdigest() + '  ' + p.relative_to(stage).as_posix() + '\n' for p in files)
(stage / 'FILES-SHA256.txt').write_text(manifest, encoding='utf-8')
if system == 'Windows':
    archive = out / (name + '.zip')
    with zipfile.ZipFile(archive, 'x', zipfile.ZIP_DEFLATED) as z:
        for p in sorted(stage.rglob('*')):
            if p.is_file():
                z.write(p, p.relative_to(out).as_posix())
else:
    archive = out / (name + '.tar.gz')
    with tarfile.open(archive, 'x:gz') as tar:
        tar.add(stage, arcname=name)
print(json.dumps({'archive': str(archive), 'sha256': hashlib.sha256(archive.read_bytes()).hexdigest(), **info}, indent=2))
