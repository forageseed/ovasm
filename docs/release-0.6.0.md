# OVASM 0.6.0

Commands are grouped into six modules, `ovasm <module> <command>`, with `-h` at every level. The
flat commands of 0.5.0 are gone; there are no aliases.

## Commands

| 0.5.0 | 0.6.0 |
|---|---|
| `run --read-type hifi` (default) | `run hifi` |
| `run --read-type sr` | `run short` |
| `recruit` | `find recruit` |
| `discover` | `find discover` |
| `panel` | `find classify` |
| `label` | `find label` |
| `correct` | `assemble correct` |
| `assemble` | `assemble graph` |
| `components` | `assemble filter` |
| `unify` | `assemble unify` |
| `evidence` | `check evidence` |
| `unzip` | `check unzip` |
| `linearize` | `solve linearize` |
| `primary` | `solve primary` |
| `chromosomes` | `solve chromosomes` |
| `pan` | `pan merge` |
| `reconstruct` | `pan reconstruct` |
| `identify` | `pan identify` |
| `sharing` | `pan sharing` |

Options are unchanged, except that `run`'s `--read-type` is now the `run short` / `run hifi`
choice. The JSON report of `pan reconstruct` says `"tool": "ovasm pan reconstruct"`.

The same inputs give byte-identical molecules.fasta, assembly.gfa and organelle.gfa with 0.5.0's
flat commands (plastid HiFi target reads; 200k short reads).

## Short-read recruitment

- `find recruit --kmer-gate`: after the first pass, the library's k-mers are counted once
  (count-min, 4 GiB) and a read extends a target only when the k-mers it would add are at 1/4 to
  10 times the target's depth. Iterative baiting from short reads had followed NUMT edges into
  the nucleus (WL354's mitochondrion grew from 0.72 M to 54 M reads; CRR1151719's to 104 M).
- `run short` recruits with `--kmer-gate --mates-together` and caps each organelle's bases
  (`--recruit-max-bases`, as for HiFi). WL354: 1.12 M mitochondrial reads, 1:37:20, 10.9 GB;
  CRR2683883: 0.61 M, 1:47:44, 9.6 GB. Both still assemble into about fifty pieces: the
  read-pair closure and frontier walk of `pan reconstruct` are not yet part of `run short`.
