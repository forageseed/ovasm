# Schoenoplectus: assembly records checked against the original reads

The Schoenoplectus mitochondrial "reference" used in earlier scoring is the Oatk data package's
`MITO.fasta`, an assembler output, so it is not independent truth. Here both record sets are
checked against the original HiFi reads with the same rules, written down before the results
were seen. Development data, already seen; not an independent biological acceptance.

## Method

- Records: OVASM `--rust-hpc --k 1001 --s 31 --min-count 80` primary output (25 records,
  2,094,440 bp; the 8 main records are compared) and Oatk `MITO.fasta` (8 linear records).
- Reads: all 1,044,826 reads of ERR8705856 (ENA), mapped separately to each set with
  `minimap2 2.30 -x map-hifi -c --cs --secondary=no`, with the same Schoenoplectus plastid sequence
  (Oatk `lpSchLacu1` plastid, 232,552 bp) added to both sets as a competing target.
- Checker: [`scripts/assembly_benchmarks/validate_records.py`](../scripts/assembly_benchmarks/validate_records.py)
  and its synthetic self-test `validate_records_selftest.py` (a correct record passes; a 10-kb
  deletion misjoin, a substitution and a chimeric join are each caught).
- A record is **read-validated** when, away from its linear ends, every position is spanned by
  at least 3 continuous alignments with 1 kb on each side (an indel of 50 bp or more breaks an
  alignment), there are no breakpoints, and no position with depth >= 5 has at least half of
  its alignments disagreeing (mismatch, deletion or insertion; each alignment counts once per
  position) or carrying a >= 50-bp indel.

Two method changes are recorded: the 50-bp indel cut was added after the synthetic self-test
showed that a misjoin can align through as one long insertion (before any real data); counting
each alignment once per position was fixed after the first real results, where a mismatch next
to an insertion was counted twice in one strand's coordinates only. Both apply to both sets.

## Results

| Main records read-validated | Mitochondrial targets only | With the plastid competing |
|---|---|---|
| OVASM (8 main records) | 5 / 8 | **8 / 8** |
| Oatk (ctg1–ctg8) | 4 / 8 | **7 / 8** |

- Continuity: no breakpoints in any main record of either set; every position is spanned by at
  least 75 reads with 1 kb on each side.
- Three high-depth regions (about 2,700–3,400x against 130–170x for the records) disagreed in
  both sets when only mitochondrial targets were offered, and all became concordant once the
  plastid was added: plastid reads on plastid-derived inserts, not assembly errors.
- Oatk ctg5 keeps disagreeing sites near 143.5 kb: about 300 reads split into two groups. The
  150 reads with the second variant all map, in the OVASM set, to a 69,079-bp sequence that
  OVASM places after the end of ctg1 and that the Oatk set does not contain.
- That 69 kb is joined to the ctg1 sequence by 93 reads spanning 1 kb on each side (46 at
  5 kb), has a median depth of 142x (record 162x, minimum 97x), at least 90 reads span every
  position, it matches the plastid in only 127 bp and carries a complete atp8 panel hit.

## What this does and does not show

On the 8 main mitochondrial records, the reads agree with OVASM everywhere and with Oatk except
where the Oatk set lacks the 69-kb sequence. It does not make OVASM's records the full truth:

- the repeat node u17 has two exits, each joined by many complete reads (140 and 125); the
  records represent one of the two configurations;
- joins inside repeats longer than the reads cannot be judged by reads;
- 17 further short OVASM records (53,864 bp, depth 500–24,000x) are mostly not validated and
  their origin is unresolved; they are not counted as mitochondrial sequence;
- physical molecule number and topology are not established.

Using these 8 records as a reference for scoring OVASM itself would be self-assessment; they can
serve as a read-validated development reference for other tools or later versions.
