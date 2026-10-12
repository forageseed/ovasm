# Short-read local adjacency evidence

OVASM previously required long-span anchors to support every assembly link. For a dense k=127 graph with a 126 bp overlap and 40 bp anchors on each side, a strict crossing needs at least 206 bp. A 150 bp read cannot meet that requirement even when it contains the exact sequence of a local edge.

Short-read mode now additionally counts the exact overlap-edge word, including a base on each side of the overlap (and at least the evidence k bases for shorter overlaps). The graph overlap must match exactly. The complete word is verified after indexed lookup; forward and reverse-complement matches count once per read. Identical words at several links, or words also present inside a graph segment, are ambiguous and assign no link.

This evidence supports local adjacency. It never enters long-span crossings, read walks, repeat traversals or phased paths. Direct edges between nonrepeat anchors can use the maximum of strict and local support, a conservative lower bound because the two read sets may overlap. Repeat entry/exit phasing retains the existing long-span requirements.

The public `run --read-type sr` workflow enables the rule; `evidence --short-read-adjacency` enables it explicitly. The Python evidence API enables it for `preset="sr"`. HiFi defaults remain unchanged.

Evidence JSON keeps `reads` and `adjacency_reads` separate, with `adjacency_ambiguous_reads` and `adjacency_supported_links`. OV-GFA keeps strict support in `ev:i`, local support in `aj:i`, and ambiguous local counts in `aa:i`. Nonzero local tags are emitted independently; an `aj` tag is not repeat phasing or physical topology evidence.

## Existing-data development check

Six fixed original graphs and complete previously recruited FASTQ sets were compared, without repeating whole-WGS recruitment. Strict link counts were unchanged in all six. In one plastid case, supported bridges increased from 60 to 603, unbridged anchors fell from 352 to 1, and candidate FASTA records fell from 37 to 3. Its 47 partial paths and `decisive=false` remain. Mitochondrial changes were small; graphs with missing edges cannot be repaired by counting evidence on existing edges.

All six graphs were separately audited for words occurring inside segments: none were found, so the additional ambiguity rule excludes no previously assigned link in this panel. These are seen-data development comparisons, not independent truth validation or proof of complete mitochondrial assemblies. Original benchmark results and input scopes are preserved.
