# Rust graph kernels

Current graph construction is implemented in Rust, with no graph-kernel FFI or C
build step. `ovasm run` and ordinary `assemble` use the base-coordinate sparse
graph; short reads use dense mode. `assemble --multichrom-graph` and `--rust-hpc`
select the experimental homopolymer-coordinate implementation. The former
defaults to k=1001 runs; the latter requires an explicit k. Both require a
read-derived `--min-count`.

HPC construction compresses homopolymers, builds the Rust syncmer graph and
restores lengths from distinct original read/run observations at shared graph
coordinates. Restored overlaps must agree exactly. Censored terminal runs do not
vote; unsupported positions fail instead of receiving invented sequence.

The HPC implementation does not yet provide read correction or repeat unzipping.
Reports mark both capabilities false. This is separate from the default Rust
assembler's existing correction, gap-rescue and original-read junction routes.
Local link support does not certify complete paths, origin or physical topology.

Syncmers are a published method ([Edgar, 2021](https://peerj.com/articles/10805/));
using them does not claim invention of syncmers or sparse de Bruijn graphs.
The implementation is in `src/assemble.rs` and `src/hpc.rs`.

```sh
cargo build --release --locked
ovasm assemble --reads sample.hifi.fastq.gz --multichrom-graph \
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

Full-input multichrom scores are recorded separately when terminal. These are
previously seen development data, not independent biological acceptance.
Version-tagged benchmark records preserve the implementation used for each run.
