# OVASM mechanism figure

![OVASM mechanism](ovasm-mechanism-v2.svg)

**Figure. Evidence-aware organelle assembly with OVASM.** Reads enter technology-specific graph construction: native closed-syncmer or portable graphs for HiFi, correction followed by graph construction for noisy long reads, and dense k-mer graphs for short reads. Seed recruitment or discovery supplies candidate organelle reads; the native multichrom HiFi strategy can construct the graph from the original whole-genome input. Read-supported flank pairings resolve repeats where evidence permits; ambiguous or partial walks remain alternatives. Native multichrom consensus pools observations at shared graph coordinates and counts each read position once. Component assessment combines marker and coverage evidence with optional original-read source contexts. Source-exclusive evidence can exclude eligible singleton cycles; compatible or unresolved components retain their evidence status. Export provides organelle sequences, graphs and evidence reports. Optional evaluation against a same-sample reference reports relaxed recovery and strict specificity, correspondence and length criteria.

Solid arrows show processing in numbered order (1 → 2 → 3 → 4 → 5 → 6). Dashed arrows show original-read evidence entering component assessment and optional evaluation after export (6 → 7). The schematic does not represent measured lengths, chromosome counts or demonstrated physical topology. Graph cycles and FASTA record counts do not establish physical molecule number or topology.

The figure describes the combined OVASM implementation. Technology selectors and evidence features depend on the execution strategy: the standalone public `run` command supports HiFi and short reads; the OrganelleVerse integration supplies the ONT/CLR adapter. Source review used public OVASM `b529d67` and OrganelleVerse integration `330ce767`. This is a mechanism diagram, not a benchmark result or a claim of independent biological validation.

## Files and reproduction

- [Editable SVG](ovasm-mechanism-v2.svg): the user-selected v2 artwork, preserved byte for byte.
- [Vector PDF](ovasm-mechanism-v2.pdf) and [600 dpi PNG](ovasm-mechanism-v2.png): matching exports.
- [Python source](draw_ovasm_mechanism_v2.py): Matplotlib source; requires Matplotlib and a compatible sans-serif font. Run `python draw_ovasm_mechanism_v2.py` to export into this directory. SVG metadata and font metrics may vary by environment.

Final dimensions are 183 × 124 mm. Server audits of the supplied PDF found a minimum text size of 5.2 pt, no text below 5 pt, and no detected text collisions or clipping. The layout is a continuous seven-step schematic rather than a quantitative multi-panel plot. The previous workflow artwork remains at `../workflow.svg`.
