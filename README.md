# ovasm

**Assemble plant mitochondrial and chloroplast genomes from sequencing reads — and get the evidence, not just one sequence.**

[English](README.md) | [简体中文](README.zh-CN.md)

Part of [OrganelleVerse](https://github.com/forageseed/organelleverse).

<p align="center"><img src="docs/workflow.svg" alt="ovasm workflow" width="780"></p>

## Why ovasm

- **One pass over the raw reads.** Organelle reads are picked out once; everything after that works on that small set.
- **An assembly graph with evidence, not just a sequence.** Every junction comes with the number of reads that support it, and `molecules.fasta` holds the representative molecule(s). When the reads cannot tell structures apart, `decisive` is `false` instead of a silent guess.
- **One small Rust program.** No conda, no containers, no system libraries (gzip is pure Rust), so it builds on Windows, Linux and macOS.
- **Pangenome-ready.** The graph is also written as OV-GFA, which `ovasm pan` merges across samples.

## Install

```bash
git clone https://github.com/forageseed/ovasm.git
cd ovasm
cargo build --release            # Rust 1.80 or newer
./target/release/ovasm --version
```

## Use

You need whole-genome reads (HiFi, or Illumina with `--read-type sr`) and a close relative's organelle genomes as seeds.
Give both organelles' seeds so each claims its own reads:

```bash
ovasm run --reads sample.hifi.fastq.gz --organelle both \
  --seeds mitochondrion=mt_reference.fasta --seeds plastid=pt_reference.fasta \
  --out out/ --threads 8
```

Results are written to `out/mitochondrion/` and `out/plastid/` (for one organelle, `--organelle mitochondrion` or
`plastid` writes straight into `out/`):

| File | Content |
|---|---|
| `assembly.gfa` | assembly graph; segments carry depth |
| `evidence.json` | reads supporting every link and branch; which segments are repeats |
| `molecules.fasta` | the representative molecule(s), unfolded along read-supported paths |
| `linearization.json` | how that representative was chosen and whether the reads can decide (`decisive`) |
| `organelle.gfa` | the graph as OV-GFA, ready for `ovasm pan` |
| `summary.json` | run summary: k used, `k_accepted`, molecule lengths, links, files |

Check `summary.json` (`k_accepted`, `decisive`) before trusting a result.

Without a reference, `--recruit discover` tries to find organelle reads by k-mer depth (HiFi only). It needs the
organelle reads to stand out from the nuclear background and can fail; give seeds when you have them.
The [OrganelleVerse](https://github.com/forageseed/organelleverse) Python package ships a ready-made seed database.

## Commands

`ovasm run` chains the stages below; each is also a command (`ovasm <command> --help`).

| Command | Does |
|---|---|
| `run` | the whole pipeline: recruit, assemble, evidence, linearize, unify |
| `recruit` | pick organelle reads out of whole-genome data using seeds |
| `discover` | find organelle reads without seeds, by k-mer depth |
| `panel` | say which organelle a set of reads comes from, by conserved proteins |
| `correct` | self-correct noisy long reads (ONT, CLR) |
| `assemble` | build the assembly graph |
| `evidence` | count read support for every link of a graph |
| `linearize` | pick the representative sequence(s) |
| `unify` | rewrite a graph as OV-GFA |
| `pan` | merge OV-GFA graphs of several samples into one pangenome graph |
| `components` | keep only the connected parts that are the target organelle |
| `identify` | find the closest reference genomes of an assembly |
| `label` | label graph pieces mitochondrion or plastid by gene content |

## Status

Tested on gold-standard genomes (Nipponbare, *Arabidopsis thaliana* Col-0, *Salvia miltiorrhiza*). Testing across
more species is still in progress, so treat results on a new species as something to check.

## Use from OrganelleVerse

```python
ov.assembly.assemble(data, organelle="mitochondrion", method="ovasm")
```

OrganelleVerse finds the program through the `ORG_VERSE_OVASM_BIN` environment variable, or on your `PATH`.

## License

MIT
