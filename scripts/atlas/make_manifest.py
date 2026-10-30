"""Build the HiFi atlas manifest: one row per species, the complete read files of one sample.

Sources, in order of preference when a species is in several (the others go to `also_in`):
  dtol      DToL HiFi (dtol_195_organelle/raw/<Species_ACCESSION>/<run>.fastq.gz)
  bigmito   large-mitochondrion HiFi (bigmito_pacbio/raw and raw_ncbi/<Species_ACCESSION>/)
  fabaceae  legume HiFi (fabaceae_pangenome/hifi and hifi_ncbi/<Species>/<biosample>/)
A read file counts only with its download marker (<file>.ok or <file>.ok.json) beside it.

usage: make_manifest.py OUT.tsv
"""
import csv
import re
import sys
from pathlib import Path

DTOL = Path('/home/jiazc/data16t-2/projects/dtol_195_organelle')
BIGMITO = Path('/home/jiazc/data16t-2/projects/bigmito_pacbio')
FABACEAE = Path('/home/jiazc/data16t-2/projects/raw_reads/fabaceae_pangenome')
READS = re.compile(r'\.(fastq|fq)(\.gz)?$')


def complete_reads(d):
    """Read files in d whose download finished: a marker <file>.ok, <file>.ok.json or
    <run>.ok.json beside it, and no aria2 control file (a download still open)."""
    out = []
    for f in sorted(d.iterdir()):
        if f.is_file() and READS.search(f.name):
            run = READS.sub('', f.name)
            done = any(Path(d, n).exists() for n in (f.name + '.ok', f.name + '.ok.json', run + '.ok.json'))
            if done and not Path(str(f) + '.aria2').exists():
                out.append(f)
    return out


def species_of(name):
    """'Carex_echinata_PRJEB88757' -> 'Carex_echinata'."""
    return re.sub(r'_(PRJ[A-Z]+\d+|SAM[A-Z]+\d+)$', '', name)


def dtol_clades():
    clades = {}
    t = DTOL / 'species_summary.tsv'
    if t.exists():
        for row in csv.DictReader(open(t), delimiter='\t'):
            clades[row['species']] = row.get('clade', '')
    return clades


def rows():
    clades = dtol_clades()
    for d in sorted((DTOL / 'raw').iterdir()):
        if d.is_dir():
            yield 'dtol', species_of(d.name), d.name, clades.get(d.name, ''), complete_reads(d)
    for sub in ('raw', 'raw_ncbi'):
        base = BIGMITO / sub
        if base.is_dir():
            for d in sorted(base.iterdir()):
                if d.is_dir():
                    yield 'bigmito', species_of(d.name), d.name, '', complete_reads(d)
    for sub in ('hifi', 'hifi_ncbi'):
        base = FABACEAE / sub
        if base.is_dir():
            for sp in sorted(p for p in base.iterdir() if p.is_dir()):
                for s in sorted(p for p in sp.iterdir() if p.is_dir()):
                    yield 'fabaceae', sp.name, f'{sp.name}_{s.name}', 'Fabaceae', complete_reads(s)


def main():
    out = Path(sys.argv[1])
    chosen, skipped = {}, []
    for source, species, sample, clade, reads in rows():
        if not reads:
            skipped.append((source, sample, 'no complete read file'))
            continue
        if species in chosen:
            chosen[species]['also_in'].append(f'{source}:{sample}')
            continue
        chosen[species] = {
            'sample': sample,
            'species': species,
            'source': source,
            'clade': clade,
            'read_type': 'hifi',
            'reads': ';'.join(str(f) for f in reads),
            'bytes': sum(f.stat().st_size for f in reads),
            'also_in': [],
        }
    fields = ['sample', 'species', 'source', 'clade', 'read_type', 'bytes', 'reads', 'also_in']
    with open(out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        w.writeheader()
        for r in sorted(chosen.values(), key=lambda r: (r['source'], r['species'])):
            w.writerow({**r, 'also_in': ','.join(r['also_in'])})
    by = {}
    for r in chosen.values():
        by[r['source']] = by.get(r['source'], 0) + 1
    print(f'{len(chosen)} species {by}; skipped {len(skipped)}', file=sys.stderr)
    for s in skipped[:20]:
        print('  skipped', *s, file=sys.stderr)


if __name__ == '__main__':
    main()
