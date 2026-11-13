# Stable whole-read protein candidate export

Use the diagnostic with a stable, per-organelle output cap.

    ovasm find protein-seeds --reads reads.fastq.gz --lineage landplants \
      --whole-read-candidates --whole-read-candidate-cap 20000 \
      --whole-read-export-order stable-record --out candidates

The diagnostic scans all records. Each organelle retains its lowest-priority record occurrences up to the requested cap. Priority hashes the full FASTX identifier, sequence and quality with the same versioned format used by stable depth-cluster sampling. Canonical FASTA record bytes break hash ties. Identical selected records appear once per selected occurrence in the FASTA. The input path, file order and record ordinal do not rank the records.

The report lists eligible, exported and cap-omitted occurrences for each organelle. Stable mode also reports distinct retained priority and crop pairs, temporary file writes and peak temporary crop payload bytes. The selector retains at most the cap in index entries for each role. Crop payloads are staged on disk and compared using fixed-size buffers. The staged payload holds at most the cap plus one incoming crop per role. The existing FASTX parser and protein scanner use their own batch memory. Temporary crop writes can increase when better-ranked records arrive late. A staging or output error leaves the diagnostic unpublished.

The diagnostic defaults to input-order export. The separate --sample-order flag controls depth-cluster sampling. Stable diagnostic candidates still need independent source evidence before downstream use. Production fallback uses its existing selection behavior.
