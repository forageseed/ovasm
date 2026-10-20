# Rust graph kernels

Current graph construction is implemented in Rust, with no graph-kernel FFI or C
build step. One engine builds the graph in one of two coordinate encodings.
`ovasm run` and ordinary `assemble` use base coordinates (sparse graph for long
reads, dense mode for short reads). `assemble --rust-hpc` selects the
experimental homopolymer-run encoding; the hidden `--multichrom-graph` is its
legacy name and does not mean a chosen chromosome count.

In the run encoding each read is homopolymer-compressed as it is read in (no
compressed-read file is written), the same backbone, edge rule, cleaning and gap
rescue build the graph, and run lengths are restored from distinct original
read/run observations at shared graph coordinates. Restored overlaps must agree
exactly; censored terminal runs do not vote; unsupported positions fail instead
of receiving invented sequence. The shared original-read junction stage then
runs on the restored base graph. `--k` counts runs (default 1001).
`--min-count` and `--min-count-fraction` act as in base coordinates.

Give a read-derived `--min-count` for whole-genome input. When it is omitted the
engine uses one tenth of the node depth (at least 3), which suits recruited
organelle reads but not whole-genome reads: on the Schoenoplectus input
(1,044,826 reads) it resolved to 3 and kept nuclear sequence, a 644-Mb run graph,
and the next stage was stopped at 271 GiB.

Read correction is not implemented in the run encoding; for HiFi the
homopolymer compression already removes the dominant error class (see the
design document). The assembly report marks read correction and read-threaded
unzipping false because `assemble` itself does neither. Repeat phases are
separate, opt-in commands that work on the restored base graph of either
encoding: `ovasm unzip` (phase table, optional `--out-gfa` path overlay and
`--materialize-gfa` split of fully resolved one-to-one repeats), `primary` /
`linearize --phases`, and `assemble --cleaning-audit-only`. Local link support
does not certify complete paths, origin or physical topology.

Syncmers are a published method ([Edgar, 2021](https://peerj.com/articles/10805/));
using them does not claim invention of syncmers or sparse de Bruijn graphs.
The implementation is in `src/assemble.rs` and `src/hpc.rs`.

```sh
cargo build --release --locked
ovasm assemble --reads sample.hifi.fastq.gz --rust-hpc \
  --min-count 80 --threads 4 --out graph.gfa --out-json graph.json
```

80 is the fixed cutoff of the named Schoenoplectus test, not a universal value.
References are used only to score the completed output.

## Development evidence

- Eight synthetic chromosomes: 8/8 distinct complete 6-kb sequences, 48,000 bp
  total, from 258 records including reverse reads, one N and two short reads.
- Existing ZM4 MT input: 29,287 reads / 516,099,590 bp, 5 nodes / 8 links;
  PT: 95,680 reads / 1,674,746,317 bp, 1 node / 1 link. Raw-read support was
  present for 8/8 MT and 1/1 PT links. Graph commands took 30.65 s / 371.85 MiB
  and 95.94 s / 576.09 MiB GNU peak RSS, respectively; recruitment is excluded.
- PT retained the 125,509-bp representative. MT exported one 307,026-bp
  representative instead of two 160,953/146,073-bp records. This does not
  establish equivalence of the complete MT path representation.
- A separate default-route comparison preserved both ZM4 GFAs byte for byte and
  all enumerated science fields against the retained release control.

- Schoenoplectus full input (1,044,826 reads, `--min-count 80`): each of the eight
  reference records is at least 99.97% covered across the candidates, with best
  single alignments at identity 1.0, but 0/8 match record-for-record. The reference is Oatk's own output, so record
  boundaries follow its unitig breaks; it is not independent truth. Details and
  the read evidence for the differences are in
  [unified-assembly-design.zh-CN.md](unified-assembly-design.zh-CN.md).

These are previously seen development data, not independent biological
acceptance.
Version-tagged benchmark records preserve the implementation used for each run.
