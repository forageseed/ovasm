# Embedded Syncasm kernel

The C and header files in this directory originate from the Syncasm graph-building
kernel from **Oatk 1.0**, Git revision
`98cd343ff7e8ff506df7fd5e6e60ceaab7657c69` (tag `1.0`), by Chenxi Zhou.
Source: https://github.com/c-zhou/oatk/tree/1.0. See `LICENSE` and file headers.

`ovasm_bridge.c` is the OrganelleVerse binding. Cargo compiles the kernel into the
OVASM Linux executable; running OVASM does not invoke an Oatk executable. The build
requires a C compiler, an archiver and zlib development headers/library. Linux is
the supported scientific provider for this kernel; other OVASM stages remain
available on other targets.

The binding enables upstream run-coordinate read correction (edit fraction 0.02),
relative arc coverage (0.35), three rounds of read-threaded graph unzipping, and
upstream base-coordinate consensus. Coverage cutoff is supplied explicitly from
the reads. The graph kernel has no access to an assembly reference, gene database,
requested chromosome count or truth scores. Organelle selection and primary
sequence representation are performed by separate OVASM stages.

Upstream source changes retain attribution and are documented here:

- `syncasm.c`, `make_ma_block`: an empty LCS contributes no blocks. Skip the
  zero-length vector append, whose `memcpy` otherwise receives null pointers.
  Full-kernel UBSan reproduced this during the eight-circle read regression.
- `levdist.c`, `wf_extend`: load eight-byte sequence blocks with `memcpy`,
  preserving the same word comparison without dereferencing misaligned
  `uint64_t` pointers. Full-kernel UBSan reproduced the invalid load with
  one ambiguous base in a long read during correction.

- `syncerr.c`, `update_syncmer_db`: original observations are counted using the
  error-correction bit of `k_mer`, instead of the strand bit of `m_pos`.
  Consensus subsequently requires an observation whose `k_mer` correction bit is
  clear. The original strand-based deletion predicate wrongly removed reverse-only
  original observations and could retain corrected-only nodes, causing a consensus
  assertion with real recruited reads at cutoff 1. The regression fixture in
  `tests/syncmer_correction_flags.c` checks original-only, corrected-only and mixed
  observations under both strand assignments.
- `syncasm.c` and `shared_consensus.h`: overlapping unitigs now restore
  homopolymer lengths on shared HPC graph coordinates. A union-find identifies
  positions already joined by the graph's exact HPC overlaps. At those positions,
  the kernel pools original observations from every covering syncmer and counts
  each `(read index, HPC read position)` once. It uses the existing rounded mean
  estimator of run length; corrected observations remain excluded. Interior
  coordinates retain the upstream per-unitig consensus. Both orientations of an
  edge receive the resulting base overlap, without clamping a wrong overlap or
  discarding an edge. The graph still has to pass OVASM's exact base-overlap
  validation. `tests/shared_run_consensus.c` covers differing independent means,
  duplicate observations, branching overlaps, corrected-observation exclusion,
  both read and graph strands, and long runs stored in the auxiliary array.
  Syncmer placements are bucketed by read, and distinct position pairs are
  accumulated with Syncasm's existing integer-key map, cleared between reads.
  This preserves the same integer sums and counts without storing and sorting
  every duplicate base observation from the entire input.

This shared-coordinate consensus is an OVASM extension, rather than an assertion
that the unmodified Oatk kernel implements the same estimator. It changes base
consensus, and therefore needs real-read sequence regression in addition to graph
invariant tests. It does not establish chromosome completeness or physical shape.

The OVASM bridge retains the kernel's existing `scg_meta_t` until it has counted
input records and summed raw lengths saved by the reader before HPC, then frees it.
The added `sr_t.original_length` includes ambiguous bases; correction does not
change it. Run-length arrays are not used to reconstruct raw input statistics.
The report therefore counts the actual assembled input, including reads too
short to produce a syncmer, without reopening and decompressing the WGS input.
