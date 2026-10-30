[English](README.md) | [简体中文](README.zh-CN.md)

<p align="center"><img src="docs/logo.png" alt="ovasm logo" width="440"></p>

# OVASM

**Rust organelle assembly: reads → mitochondrial and plastid sequences, GFA and read evidence.**

OVASM recruits organelle reads, builds graphs and retains read-supported alternative connections. HiFi and Illumina inputs are supported; species reference graphs can assist WGS reconstruction.

![OVASM mechanism](docs/figures/ovasm-mechanism-v2.svg)

## Install

**Precompiled 0.6.0 — no compiler required.**

- **Windows x86_64:** [download ZIP](https://github.com/forageseed/ovasm/releases/download/v0.6.0/ovasm-v0.6.0-windows-x86_64.zip), extract and double-click `install.cmd`. Reopen the terminal.
- **Linux x86_64, glibc 2.34+:**

```sh
curl -fsSL https://raw.githubusercontent.com/forageseed/ovasm/v0.6.0/scripts/install.sh | sh
export PATH="$HOME/.local/bin:$PATH"
ovasm --version
```

[Downloads](https://github.com/forageseed/ovasm/releases/tag/v0.6.0) · [Offline installation](docs/installation.md)

## Use

```sh
ovasm run hifi --reads sample.hifi.fastq.gz --organelle both --threads 8 --out result
```

Use `ovasm run short` for Illumina. The output contains representative FASTA, GFA, read-support reports and source candidates. Inspect unresolved paths before interpreting graph cycles as physical molecules.

Commands are `ovasm <module> <command>`; `-h` works at every level:

| Module | Commands |
|---|---|
| `run` | `short`, `hifi`: one step from reads to graph and molecules |
| `find` | `recruit`, `discover`, `classify`, `label`: the organelle's reads |
| `assemble` | `correct`, `graph`, `filter`, `unify`: reads to an assembly graph |
| `check` | `evidence`, `unzip`: the graph against the reads |
| `solve` | `linearize`, `primary`, `chromosomes`: sequences and chromosomes |
| `pan` | `merge`, `reconstruct`, `identify`, `sharing`: many samples |

[Species graph → individual WGS](docs/reference-reconstruction.md) · [Commands](https://github.com/forageseed/ovasm/tree/main/docs)

## Results

| Complete ZM4 HiFi input, 8 threads | Wall time | Peak sampled process-group RSS |
|---|---:|---:|
| OVASM | 210.20 min | 1.983 GiB |
| OATK, including annotation | 942.35 min | 310.416 GiB |

OVASM produced MT representatives of 160,953 and 146,073 bp and a PT representative of 125,509 bp. These are sequence representatives, not certified physical chromosome counts. The table records the named full-input benchmark versions; it is not a new 0.3.0 speed measurement.

The 20-WGS first-attempt panel completed 14 commands; six hit resource limits. Complete MT assembly was not established by those execution results.

[Benchmark details](docs/benchmarks.md) · [WGS panel](docs/existing-wgs-benchmark.md) · [GFA/ONT evidence](docs/hifi-configuration-comparison.md) · [Rust HPC validation](docs/rust-kernel.md)

## Contributors and license

Codex (OpenAI): implementation, testing, figures and documentation. **MIT**.
