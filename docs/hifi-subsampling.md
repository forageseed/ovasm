# Deterministic nested HiFi subsets

Run from the formal repository, choosing a new output directory on the approved data24 volume:

```sh
python3 -I scripts/subsample_hifi.py \
  --reads /data/sample/runA.fastq.gz --source-id runA \
  --reads /data/sample/runB.fastq.gz --source-id runB \
  --seed 17 --seed 29 --seed 43 \
  --fraction 0.05 --fraction 0.10 --fraction 0.25 \
  --out /home/jiazc/data24t/NEW_APPROVED_BENCHMARK/subsets
```

The paths above are examples; the tool does not choose, download, or discover biological inputs. Give all assemblers the **same generated subset files**, rather than asking each tool to resample independently.

The tool scans the source files once, decompressing each once. Each full read gets a uniform 64-bit BLAKE2b value based on seed, source ID, one-based source record index and read ID. Thresholds select nested subsets within each seed. Different seeds are independent sampling replicates; they are not nested with one another. Duplicate read IDs are disambiguated by source ID and record index. Preserve record ordering and explicit source IDs to reproduce the same subset after moving inputs. Defaults are input ordinal plus basename; use explicit IDs for publication runs.

Fractions are inclusion probabilities, **not exact numbers of bases or measured nuclear coverage**. Whole reads are retained without truncation to a target base count. Record actual reads, bases and any independently estimated nuclear depth in the benchmark table.

Sequence and quality contents are preserved, including wrapped FASTQ. A missing final newline is added only when necessary to concatenate records from multiple sources safely. Gzip metadata is deterministic (no filename or timestamp). Compressed bytes are reproducible in the same Python/zlib environment; the uncompressed SHA-256 is the portable identity.

The output directory must not exist. While running, or if a failure occurs, it contains `INCOMPLETE`; do not benchmark such outputs. Successful completion removes that marker and writes `status: complete` in `manifest.json`. On failure, partial output is retained for diagnosis and never represented as complete. Use a new directory for a rerun.

The manifest records parameters, exact key construction, source paths/IDs/size/mtime, SHA-256 of uncompressed source records, output SHA-256, actual base fractions, counts, min/max/mean read lengths and full length histograms. Input size/mtime changes during a run are rejected. Source checksum refers to uncompressed FASTQ bytes, not compressed archive bytes.

Disk usage grows with the sum of requested fractions across seeds. All outputs remain open during the scan; keep the number of subsets modest. Compression and hashing are CPU work and should use the experiment's agreed CPU/nice limits externally.

Validation: `python3 -I scripts/test_subsample_hifi.py` passed four fixture tests covering nesting, repeated IDs, multi-source gz/plain and wrapped FASTQ, stable output across repeated runs and relocated paths, different seeds, complete read/quality lengths, overwrite rejection, incomplete failure marking and missing-final-newline handling. No real input library was scanned during implementation.

Disk guard: `--min-free-gib` defaults to 30 GiB. Free output-disk space is checked at start, after each 64 MiB of uncompressed records or 1000 records (whichever comes first), and before completion. Falling below the floor aborts and preserves `INCOMPLETE`. The floor is a periodic guard, not a space reservation; concurrent writers can still consume space between checks. Six fixture tests include startup and midstream disk exhaustion.

Compression: `--compresslevel` accepts 0–9 and defaults to 1 for throughput. The level is recorded in the manifest and affects compressed size/bytes only; selected reads and uncompressed checksums are unchanged. A seventh fixture test compares levels 1 and 9.
