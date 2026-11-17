//! `ovasm run ont` and `ovasm run hybrid`: noisy long reads (ONT, PacBio CLR), alone or with the
//! sample's short reads, in to organelle graph and molecules out.
//!
//! Noisy reads cannot extend a bait by themselves: their error k-mers are short random words that
//! the nuclear genome holds too, so `find recruit --preset ont` takes one pass and never extends
//! (see its `--iterations`). The flow extends from corrected reads instead:
//!
//! 1. one recruitment pass with the seeds (`find recruit --preset ont`, k 21: at k 15 the built-in
//!    seed library's k-mers hit every read by chance, see `RECRUIT_K`);
//! 2. self-correction of each organelle's reads (`assemble correct`, k 17 then 25);
//! 3. further rounds (`--recruit-rounds`): the round's reads and its partly covered ones (a 2 kb
//!    block of the bait, see `split_partial`), corrected together, join the seeds as bait
//!    (`find recruit --seed-reads`), so the reads that reach past the seeds, the sample's own
//!    intergenic sequence, are recruited; a round that adds less than `ROUND_GROWTH` of the
//!    bases is not used and ends the rounds. The last round's reads, corrected on their own
//!    (without the partly covered reads, which hold the NUMTs), are the ones assembled;
//! 4. `run hybrid` only: the corrected long reads seed one pass over the short reads (`find
//!    recruit --preset sr --seed-reads`, mates together), and the last round's long reads are
//!    corrected again with those short reads (`assemble correct --short-reads`, then self rounds
//!    for what they do not cover). The short reads cannot come first: without the sample's own
//!    long-read sequence as bait, seeds recruit short reads of the conserved genes only;
//! 5. the corrected reads go through the HiFi stages as target reads (`run hifi --read-set
//!    target`: sparse graph at the automatic k with the k ladder, evidence, unringing, unified
//!    graph, chromosomes, `check assess`);
//! 6. `run hybrid` only: the short reads check the same graph as well (`evidence.short.json`,
//!    short-read adjacency), beside the corrected long reads' `evidence.json`.
//!
//! Every stage leaves its usual report; `summary.json` is that of `run hifi` with
//! `read_type: ont | ont_hybrid` and a `long_read` record of the rounds, corrections and short
//! reads.

use std::collections::{BTreeMap, BTreeSet};
use std::ffi::OsString;
use std::fs;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{bail, Context, Result};
use serde_json::{json, Map, Value};

use crate::pipeline::{
    self, evidence_argv, read_json, write_json, Argv, Mode, OrganelleChoice, ReadSet, ReadType,
    MIN_RECRUITED_BASES,
};
use crate::bait::BaitSet;
use crate::kmer::CanonicalKmers;
use crate::recruit::{classify, peak_rss_mb, Call, ClassifyParams};

/// k of long-read recruitment. `find recruit --preset ont` defaults to 15, which with the
/// built-in seed library (tens of Mb of k-mers, a twentieth of all 15-mers) covers every read
/// by chance: ZM4's ONT library (CRR2374622, R10.4.1) gave 61 Gb of "plastid" reads out of 131 Gb.
/// At 21 a chance hit is a hundred thousand times rarer, and a read of 5% error still keeps a
/// third of its 21-mers.
pub const RECRUIT_K: usize = 21;

/// The k of corrected noisy reads' graph is `ASSEMBLE_K_HIGH` unless the graph at `ASSEMBLE_K_LOW`
/// holds `ASSEMBLE_K_GROWTH` times its sequence or more, which means the higher k lost sequence to
/// errors the correction left (Nipponbare ONT: 273 kb at 201, 405 kb at 101, of 376 kb; ZM4 ONT:
/// 299.7 kb at 201, 300.2 kb at 101).
pub const ASSEMBLE_K_R10: usize = 501;
pub const ASSEMBLE_K_HIGH: usize = 201;
pub const ASSEMBLE_K_LOW: usize = 101;
pub const ASSEMBLE_K_GROWTH: f64 = 1.10;

fn parse_assemble_k(s: &str) -> Result<usize, String> {
    let k: usize = s.parse().map_err(|_| format!("{s:?} is not a number"))?;
    if crate::assemble::K_LADDER.contains(&k) {
        Ok(k)
    } else {
        Err(format!("k must be one of {:?}", crate::assemble::K_LADDER))
    }
}

/// A recruitment round that adds less than this share of the previous round's bases is not
/// used, and ends the rounds.
pub const ROUND_GROWTH: f64 = 0.05;

#[derive(clap::Args)]
pub struct LongOptions {
    /// Organelle to assemble; `both` recruits once and assembles each into <out>/<organelle>/.
    #[arg(long, value_enum)]
    organelle: OrganelleChoice,
    /// Output directory: what `run hifi` writes (assembly.gfa, evidence.json, molecules.fasta,
    /// organelle.gfa, summary.json ...), plus recruit/round<N>/ and correct/round<N>/ (and
    /// recruit/short/, correct/hybrid/ for `run hybrid`).
    #[arg(long)]
    out: PathBuf,
    /// whole-genome: recruit the organelle's reads first. target: the long reads are the
    /// organelle's already and are only corrected (one organelle; `run hybrid` still recruits
    /// the short reads).
    #[arg(long, value_enum, default_value_t = ReadSet::WholeGenome)]
    read_set: ReadSet,
    /// Seed FASTA as ORGANELLE=PATH (repeat for more files); a bare PATH seeds the --organelle.
    /// Without --seeds the built-in seed library of --lineage is used, written to <out>/seeds/.
    /// Both organelles are seeded unless --single-seed, so that each claims its own reads.
    #[arg(long = "seeds")]
    seeds: Vec<String>,
    /// Built-in seed library: landplants (mitochondrion and plastid), metazoa or fungi
    /// (mitochondrion only).
    #[arg(short = 'x', long, value_enum, default_value_t = crate::seeddb::Lineage::Landplants)]
    lineage: crate::seeddb::Lineage,
    /// Recruit one organelle without seeding the other (not advised).
    #[arg(long)]
    single_seed: bool,
    /// Recruitment rounds over the long reads: 1 recruits with the seeds only; each further
    /// round adds the previous round's corrected reads to the bait (each round reads the whole
    /// input again). Rounds stop early once one adds less than a twentieth of the bases.
    #[arg(long, default_value_t = 2)]
    recruit_rounds: usize,
    /// k of long-read recruitment (`find recruit --k`); see `longread::RECRUIT_K`.
    #[arg(long, default_value_t = RECRUIT_K)]
    recruit_k: usize,
    /// Keep at most about this many recruited long-read bases per organelle, which bounds the
    /// time and memory of correction (it holds the reads in memory); 0 for no cap. 400 Mb is a
    /// 300 kb mitochondrion at 1,300x.
    #[arg(long, default_value_t = 400_000_000)]
    recruit_max_bases: u64,
    /// k the graph of the corrected reads starts at (odd, one of 1001, 701, 501, 351, 251, 201,
    /// 151, 101); it steps down as `run hifi` does when the molecules explain under half of the
    /// graph. Default: 201, or 101 when the graph at 101 holds a tenth more sequence than at 201
    /// (`choose_assemble_k`). Corrected noisy reads keep errors that repeat across reads (ONT
    /// homopolymers), and at the HiFi default of 1001 their graph breaks (Nipponbare ONT, 600x:
    /// 1 of 16 junctions at 1001, 9 at 201, 11 at 101).
    #[arg(long, value_parser = parse_assemble_k)]
    assemble_k: Option<usize>,
    /// fast or standard, as for `run hifi` (standard: downsampled replicates of the corrected
    /// reads).
    #[arg(long, value_enum, default_value_t = Mode::Fast)]
    mode: Mode,
    /// Salt of every random draw.
    #[arg(long, default_value_t = 0)]
    salt: u64,
    /// Sample name in the unified graph (segment names, PanSN paths).
    #[arg(long, default_value = "sample")]
    sample: String,
    /// Profile organelle base modifications (`CG`/`CHG`/`CHH` 5mC and `6mA`) and NUMT/NUPT separation from modBAM inputs.
    #[arg(long)]
    methylation: bool,
    /// Explicit modBAM file(s) with `MM`/`ML` tags for `--methylation` when `--reads` are FASTQ/FASTA.
    #[arg(long = "methylation-bam")]
    methylation_bam: Vec<PathBuf>,
    /// Maximum per-read CG 5mC fraction allowed during BAM read ingestion (e.g. 0.40 to drop hypermethylated NUMTs/NUPTs).
    #[arg(long)]
    max_read_methylation: Option<f64>,
    #[arg(short = 't', long, value_parser = clap::value_parser!(u16).range(1..=256))]
    threads: Option<u16>,
    /// Set by `run_ont` / `run_hybrid`: `run ont` defaults to `Auto` (`R10` phased overlap EC on
    /// R10 FASTQ, `Legacy` otherwise); `run hybrid` always uses `Legacy`.
    #[arg(skip = crate::r10_correct::OntChemistry::Auto)]
    chemistry: crate::r10_correct::OntChemistry,
}

#[derive(clap::Args)]
pub struct OntArgs {
    /// Long read files (ONT or PacBio CLR; FASTQ/FASTA, optionally .gz); repeat for several.
    #[arg(long = "reads", required = true)]
    reads: Vec<PathBuf>,
    /// ONT error-correction chemistry: `auto` selects hifiasm-style phased overlap correction on
    /// R10 FASTQ reads and preserves the `k = 17, 25` solid-kmer dBG scheme for R9, CLR, and FASTA.
    #[arg(
        long,
        value_enum,
        default_value_t = crate::r10_correct::OntChemistry::Auto,
        conflicts_with = "r10"
    )]
    chemistry: crate::r10_correct::OntChemistry,
    /// Shorthand for `--chemistry r10` (ONT R10.4 / Simplex phased overlap error correction).
    #[arg(long)]
    r10: bool,
    #[command(flatten)]
    options: LongOptions,
}

#[derive(clap::Args)]
pub struct HybridArgs {
    /// Long read files of the sample (ONT or PacBio CLR); repeat for several.
    #[arg(long = "long", required = true)]
    long: Vec<PathBuf>,
    /// Short read files of the same sample (Illumina, whole-genome or organelle reads); repeat
    /// for several, e.g. both mates.
    #[arg(long = "short", required = true)]
    short: Vec<PathBuf>,
    /// Keep at most about this many recruited short-read bases per organelle (they are held in
    /// memory for the hybrid correction); 0 for no cap.
    #[arg(long, default_value_t = 1_000_000_000)]
    short_max_bases: u64,
    #[command(flatten)]
    options: LongOptions,
}

/// The short reads of `run hybrid`.
struct Short<'a> {
    reads: &'a [PathBuf],
    max_bases: u64,
}

fn log(message: &str) {
    eprintln!("[ovasm] run: {message}");
}

pub fn run_ont(mut args: OntArgs) -> Result<()> {
    args.options.chemistry = crate::r10_correct::OntChemistry::resolve(args.chemistry, args.r10);
    finish(run(&args.reads, None, &args.options)?, &args.options)
}

pub fn run_hybrid(mut args: HybridArgs) -> Result<()> {
    args.options.chemistry = crate::r10_correct::OntChemistry::Legacy;
    let short = Short {
        reads: &args.short,
        max_bases: args.short_max_bases,
    };
    finish(run(&args.long, Some(short), &args.options)?, &args.options)
}

fn finish(report: Value, o: &LongOptions) -> Result<()> {
    eprintln!(
        "[ovasm] run: summary -> {}; peak RSS {}",
        o.out.join("summary.json").display(),
        report["peak_rss_mb"]
            .as_f64()
            .map_or("n/a".into(), |m| format!("{m:.0} MB"))
    );
    Ok(())
}

/// ORGANELLE=PATH (or a bare PATH for the single --organelle) into seed files per organelle.
fn seed_map(specs: &[String], single: Option<&str>) -> Result<BTreeMap<String, Vec<PathBuf>>> {
    let mut seeds: BTreeMap<String, Vec<PathBuf>> = BTreeMap::new();
    for spec in specs {
        let (name, path) = match (spec.split_once('='), single) {
            (Some((name, path)), _) => (name, path),
            (None, Some(organelle)) => (organelle, spec.as_str()),
            (None, None) => {
                bail!("--seeds must be ORGANELLE=PATH with --organelle both, got {spec:?}")
            }
        };
        if name != "mitochondrion" && name != "plastid" {
            bail!("--seeds: organelle must be mitochondrion or plastid, got {name:?}");
        }
        if !Path::new(path).is_file() {
            bail!("seed FASTA not found: {path}");
        }
        seeds.entry(name.to_string()).or_default().push(path.into());
    }
    Ok(seeds)
}

fn validate(long: &[PathBuf], short: Option<&Short>, o: &LongOptions) -> Result<()> {
    for r in long.iter().chain(short.map_or(&[][..], |s| s.reads)) {
        if !r.is_file() {
            bail!("reads not found: {} (a FASTQ/FASTA file, optionally .gz)", r.display());
        }
    }
    if o.recruit_rounds == 0 {
        bail!("--recruit-rounds must be at least 1");
    }
    if !(13..=crate::kmer::MAX_K).contains(&o.recruit_k) {
        bail!("--recruit-k must be between 13 and {}", crate::kmer::MAX_K);
    }
    if o.read_set == ReadSet::Target && o.organelle == OrganelleChoice::Both {
        bail!(
            "--read-set target takes the long reads of one organelle; run mitochondrion and \
             plastid each on their own reads"
        );
    }
    Ok(())
}

/// A partly covered read joins its target's bait only when its depth is at least this share of
/// the depth of the reads the target recruited. ZM4's ONT: the plastid seed is close, so the
/// mitochondrial reads holding plastid-derived (MTPT) stretches fitted the plastid better than
/// the distant mitochondrial seeds, and as the plastid's partly covered reads (976, 118 Mb) made
/// its bait 1.3 M k-mers for a 125 kb genome; those reads are at the mitochondrion's depth, a
/// twentieth of the plastid's. A read that reaches past its own target's bait is at its depth.
pub const PARTIAL_DEPTH_SHARE: f64 = 0.25;

/// k of the depth comparison: short enough for a read with 5-10% errors to keep many k-mers.
const DEPTH_K: usize = 17;

/// A target that recruited nothing has no depth to compare with; its partly covered reads are
/// compared with one another instead. Its own reads overlap one another (they hold a few kb of
/// one small genome), a NUMT read's nuclear flank is single-copy among them, so the reference is
/// this quantile of their depths, the organelle's reads. ZM4's ONT, 38 Gb: 481 ungated partly
/// covered mitochondrial reads as bait made 199,750 of the 400,000 reads of round 2 partly
/// covered (21 Gb, two hours of `find recruit`'s source candidates).
pub const SELF_DEPTH_QUANTILE: f64 = 0.9;

/// Reads go into `assemble correct` down to this depth (median count of their `DEPTH_K`-mers,
/// below read depth by the errors). Deeper, systematic errors recur in enough reads to pass the
/// solid count: ZM4's plastid at about 3,000x (400 Mb, the `--recruit-max-bases` cap) had a
/// solid peak of 17 at k 25 against 1,074 at 1,000x, took 2.4 hours to correct, and as bait
/// held 370,000 k-mers for a 125 kb genome.
pub const CORRECT_MAX_DEPTH: f64 = 300.0;

fn read_seqs(path: &Path) -> Result<Vec<Vec<u8>>> {
    let mut out = Vec::new();
    let meta = fs::metadata(path).with_context(|| format!("cannot access {}", path.display()))?;
    anyhow::ensure!(meta.is_file(), "not a regular file: {}", path.display());
    if meta.len() == 0 { return Ok(out); }
    let mut reader = needletail::parse_fastx_file(path)
        .with_context(|| format!("cannot read {}", path.display()))?;
    while let Some(record) = reader.next() {
        out.push(record.with_context(|| format!("malformed record in {}", path.display()))?.seq().into_owned());
    }
    Ok(out)
}

/// An eighth of the canonical `DEPTH_K`-mers of `seq`, chosen by hash (the same in every read).
fn depth_kmers(seq: &[u8]) -> impl Iterator<Item = u64> + '_ {
    CanonicalKmers::new(seq, DEPTH_K)
        .map(|(_, kmer)| kmer)
        .filter(|kmer| kmer.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 61 == 0)
}

/// Per read, the median count of its `depth_kmers` among all of `seqs` (`None` for a read too
/// short to hold one).
fn read_depths(seqs: &[&[u8]]) -> Vec<Option<u32>> {
    use rayon::prelude::*;
    let mut counts: rustc_hash::FxHashMap<u64, u32> = Default::default();
    for seq in seqs {
        for kmer in depth_kmers(seq) {
            *counts.entry(kmer).or_insert(0) += 1;
        }
    }
    seqs.par_iter()
        .map(|seq| {
            let mut c: Vec<u32> = depth_kmers(seq).map(|kmer| counts[&kmer]).collect();
            if c.is_empty() {
                return None;
            }
            let mid = c.len() / 2;
            Some(*c.select_nth_unstable(mid).1)
        })
        .collect()
}

/// The `quantile` (0 to 1) of the known depths.
fn depth_quantile(depths: &[Option<u32>], quantile: f64) -> Option<f64> {
    let mut d: Vec<u32> = depths.iter().flatten().copied().collect();
    if d.is_empty() {
        return None;
    }
    let i = ((d.len() - 1) as f64 * quantile) as usize;
    Some(*d.select_nth_unstable(i).1 as f64)
}

/// Which of the `partial` reads are at the depth of the `recruited` ones: a read's depth is the
/// median count of its k-mers among both sets, and it is kept when that is at least
/// `PARTIAL_DEPTH_SHARE` of the recruited reads' median depth. Nothing recruited (a seed that
/// covers no read by 80%): of the `SELF_DEPTH_QUANTILE` of the partly covered reads' own depths.
pub fn at_target_depth(recruited: &[Vec<u8>], partial: &[&[u8]]) -> Vec<bool> {
    let all: Vec<&[u8]> =
        recruited.iter().map(Vec::as_slice).chain(partial.iter().copied()).collect();
    let depths = read_depths(&all);
    let (own, theirs) = depths.split_at(recruited.len());
    let target = depth_quantile(own, 0.5).or_else(|| depth_quantile(theirs, SELF_DEPTH_QUANTILE));
    let Some(target) = target else {
        return vec![true; partial.len()];
    };
    theirs
        .iter()
        .map(|d| d.is_some_and(|d| d as f64 >= PARTIAL_DEPTH_SHARE * target))
        .collect()
}

/// Filter partly covered reads assigned to target `t` against competing organelles that are >2x
/// deeper after normalizing each target's k-mer counts by `1 / sampling_fraction`.
///
/// When `plastid` is capped by `--recruit-max-bases` (e.g. `sampling_fraction = 0.216` on
/// Nipponbare ONT), raw counts in `plastid.fastq` understate true plastid depth nearly fivefold,
/// while noisy pure-plastid reads with an error cluster can fall just under `min_cover` on the
/// plastid seed and match a mitochondrial MTPT block (`Call::NumtCandidate(mitochondrion)`).
/// Conversely, a genuine mitochondrial read crossing an MTPT boundary (`2.5 kb mt_own + 6 kb MTPT`)
/// also has its whole-read median k-mer inside the MTPT at plastid depth, so a whole-read median
/// comparison alone would discard real crossing reads.
///
/// A read whose normalized median depth in a deeper target `u` exceeds `2 * max(med_t, ref_t)` is
/// therefore kept iff it carries a gap-bridged `>= min_own_block` span of `depth_kmers` (`>= 8`
/// sampled k-mers) that stand at target `t`'s own depth (`>= PARTIAL_DEPTH_SHARE * ref_t`) and
/// are NOT shared with any >2x deeper target `u` (`norm_u < PARTIAL_DEPTH_SHARE * ref_t` and
/// `2 * norm_u < norm_t`). Pure donor leaks have zero target-unique blocks; genuine crossing reads
/// carry their target-unique flank and are preserved.
pub fn filter_cross_organelle_partials(
    recruited: &[Vec<Vec<u8>>],
    sampling_fractions: &[f64],
    partial: &[(&[u8], usize)],
    min_own_block: usize,
) -> Vec<bool> {
    use rayon::prelude::*;
    let n = recruited.len();
    if n <= 1 || partial.is_empty() {
        return vec![true; partial.len()];
    }
    let scales: Vec<f64> = (0..n)
        .map(|t| {
            let f = sampling_fractions.get(t).copied().unwrap_or(1.0).clamp(1e-6, 1.0);
            1.0 / f
        })
        .collect();
    let mut counts: Vec<rustc_hash::FxHashMap<u64, u32>> = vec![Default::default(); n];
    let mut ref_depth: Vec<Option<f64>> = vec![None; n];

    // 1. Build clean count maps from recruited reads for every target that recruited reads.
    for t in 0..n {
        if recruited[t].is_empty() {
            continue;
        }
        for seq in &recruited[t] {
            for km in depth_kmers(seq) {
                *counts[t].entry(km).or_insert(0) += 1;
            }
        }
        let own_meds: Vec<Option<u32>> = recruited[t]
            .par_iter()
            .map(|seq| {
                let mut c: Vec<u32> = depth_kmers(seq)
                    .filter_map(|km| counts[t].get(&km).copied().filter(|&v| v >= 2))
                    .collect();
                if c.is_empty() {
                    return None;
                }
                let mid = c.len() / 2;
                Some(*c.select_nth_unstable(mid).1)
            })
            .collect();
        ref_depth[t] = depth_quantile(&own_meds, 0.5).map(|d| d * scales[t]);
    }

    // 2. Fallback for a target that recruited 0 full reads in round 1 (distant multi-species seed):
    // count its partial candidates over k-mers absent from already-recruited competing organelles.
    for t in 0..n {
        if !recruited[t].is_empty() {
            continue;
        }
        let mine: Vec<&[u8]> = partial
            .iter()
            .filter_map(|&(seq, u)| (u == t).then_some(seq))
            .collect();
        if mine.is_empty() {
            continue;
        }
        for seq in &mine {
            for km in depth_kmers(seq) {
                let in_other = (0..n).any(|u| {
                    u != t
                        && !recruited[u].is_empty()
                        && counts[u].get(&km).copied().unwrap_or(0) > 1
                });
                if !in_other {
                    *counts[t].entry(km).or_insert(0) += 1;
                }
            }
        }
        let meds: Vec<Option<u32>> = mine
            .par_iter()
            .map(|seq| {
                let mut c: Vec<u32> = depth_kmers(seq)
                    .filter_map(|km| counts[t].get(&km).copied().filter(|&v| v >= 2))
                    .collect();
                if c.is_empty() {
                    return None;
                }
                let mid = c.len() / 2;
                Some(*c.select_nth_unstable(mid).1)
            })
            .collect();
        ref_depth[t] = depth_quantile(&meds, SELF_DEPTH_QUANTILE).filter(|&d| d > 0.0);
    }

    let gap_limit = (min_own_block / 3).max(400);
    partial
        .par_iter()
        .map(|&(seq, t)| {
            if t >= n {
                return false;
            }
            let Some(rt) = ref_depth[t].filter(|&d| d > 0.0) else {
                return true;
            };
            // Measure copy depth over matched solid k-mers (`v >= 2`) rather than 0-count error
            // k-mers: once ONT read error exceeds ~4.1%, `(1 - e)^17 < 0.5` collapses an
            // unfiltered whole-read median to 0, whereas the median over matched solid k-mers
            // stays invariant to sequencing error rate and chimeric nuclear tails.
            let solid_med = |u: usize, min_hits: usize| -> Option<f64> {
                let mut c: Vec<u32> = depth_kmers(seq)
                    .filter_map(|km| counts[u].get(&km).copied().filter(|&v| v >= 2))
                    .collect();
                if c.len() < min_hits {
                    return None;
                }
                let mid = c.len() / 2;
                Some(*c.select_nth_unstable(mid).1 as f64 * scales[u])
            };
            let med_t = solid_med(t, 4).unwrap_or(0.0);
            let dominated = (0..n).any(|u| {
                u != t
                    && ref_depth[u].is_some_and(|ru| ru > 2.0 * rt)
                    && solid_med(u, 8).is_some_and(|mu| mu > 2.0 * med_t.max(rt))
            });
            if !dominated {
                return true;
            }
            // Rescue genuine crossing reads that carry a >= min_own_block target-unique flank.
            let (mut start, mut last_end, mut hits) = (0usize, 0usize, 0usize);
            for (pos, km) in CanonicalKmers::new(seq, DEPTH_K) {
                if km.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 61 != 0 {
                    continue;
                }
                let raw_t = counts[t].get(&km).copied().unwrap_or(0);
                let nt = raw_t as f64 * scales[t];
                let donor_shared = (0..n).any(|u| {
                    u != t
                        && ref_depth[u].is_some_and(|ru| ru > 2.0 * rt)
                        && {
                            let nu = counts[u].get(&km).copied().unwrap_or(0) as f64 * scales[u];
                            nu > 2.0 * rt && nu > 2.0 * nt
                        }
                });
                if donor_shared {
                    // Donor-shared stretches break target-unique runs so scattered SNPs inside a
                    // diverged MTPT cannot chain across conserved donor sequence (e.g. NUMT+MTPT).
                    hits = 0;
                    continue;
                }
                let own_only = raw_t >= 2
                    && nt >= PARTIAL_DEPTH_SHARE * rt
                    && (0..n).all(|u| {
                        if u == t || !ref_depth[u].is_some_and(|ru| ru > 2.0 * rt) {
                            return true;
                        }
                        let nu = counts[u].get(&km).copied().unwrap_or(0) as f64 * scales[u];
                        nu < PARTIAL_DEPTH_SHARE * rt && nu * 2.0 < nt
                    });
                if !own_only {
                    continue;
                }
                if hits == 0 || pos > last_end + gap_limit {
                    start = pos;
                    hits = 1;
                } else {
                    hits += 1;
                }
                last_end = pos + DEPTH_K;
                if last_end - start >= min_own_block && hits >= 8 {
                    return true;
                }
            }
            false
        })
        .collect()
}

/// Short tandem repeats masked out of corrected reads used as bait: a run of at least this many
/// bases of period 1 to `TANDEM_MAX_PERIOD`. Microsatellites occur all over the nuclear genome,
/// so their k-mers in a bait recruit nuclear reads that are mostly microsatellite: ZM4's ONT at
/// 12 Gb recruited two such reads into the mitochondrion's 81 in round 2, and after hybrid
/// correction their (GA)n set the depth of the automatic k (23,504 at k 1001), which left the
/// graph without a segment. The run is one recruitment k-mer (`RECRUIT_K`): correction leaves
/// errors in a long microsatellite, which cut it into runs of 20 to 30 bases, and each of those
/// still held a k-mer of pure repeat.
pub const TANDEM_MIN_RUN: usize = RECRUIT_K;
const TANDEM_MAX_PERIOD: usize = 6;

/// Where `seq` repeats itself at a period of 1 to `TANDEM_MAX_PERIOD` for `TANDEM_MIN_RUN` bases or
/// more, as N.
pub fn mask_tandem(seq: &mut [u8]) -> usize {
    let mut mask = vec![false; seq.len()];
    for p in 1..=TANDEM_MAX_PERIOD {
        let mut start = 0;
        for i in 0..=seq.len().saturating_sub(p) {
            let same = i + p < seq.len()
                && seq[i].eq_ignore_ascii_case(&seq[i + p])
                && seq[i] != b'N';
            if !same {
                // seq[start..i + p] repeats at period p
                if i + p - start >= TANDEM_MIN_RUN && i > start {
                    mask[start..i + p].iter_mut().for_each(|m| *m = true);
                }
                start = i + 1;
            }
        }
    }
    let mut n = 0;
    for (b, m) in seq.iter_mut().zip(mask) {
        if m {
            *b = b'N';
            n += 1;
        }
    }
    n
}

/// `src` (FASTA or FASTQ) to the FASTA `dst` with its short tandem repeats masked
/// (`mask_tandem`); the bases masked.
fn mask_tandem_fasta(src: &Path, dst: &Path) -> Result<u64> {
    let mut w = BufWriter::new(fs::File::create(dst)?);
    let mut masked = 0u64;
    if src.is_file() && fs::metadata(src)?.len() > 0 {
        let mut reader = needletail::parse_fastx_file(src)
            .with_context(|| format!("cannot read {}", src.display()))?;
        while let Some(record) = reader.next() {
            let record = record.with_context(|| format!("malformed record in {}", src.display()))?;
            let mut seq = record.seq().into_owned();
            masked += mask_tandem(&mut seq) as u64;
            w.write_all(b">")?;
            w.write_all(record.id())?;
            w.write_all(b"\n")?;
            w.write_all(&seq)?;
            w.write_all(b"\n")?;
        }
    }
    w.flush()?;
    Ok(masked)
}

/// A number in [0, 1) for a read name and a salt, uniform over the names of real runs: FNV-1a over
/// the name, then the splitmix64 finaliser. Samples the reads a share is kept of (`cap_depth`,
/// the partial reads' cap). `recruit::unit_hash` does not do this for Nanopore names
/// (`SRR25241091.96/1 <uuid> length=37207`): over the 32,383 reads of Nipponbare ONT's plastid
/// its values lay in [0, 0.3), so keeping a share of 0.16 kept 75% of the reads, the plastid
/// went into correction at 2,240x and its graph came out empty.
fn read_share(id: &[u8], salt: u64) -> f64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in id {
        h = (h ^ u64::from(b)).wrapping_mul(0x0000_0100_0000_01b3);
    }
    let mut z = h ^ salt.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^= z >> 31;
    (z >> 11) as f64 / (1u64 << 53) as f64
}

/// `reads` down to `CORRECT_MAX_DEPTH`: when their median depth (see `read_depths`) is above it,
/// the reads whose id hashes below the share are written to `out` and it is returned in their
/// place. Also the record of it: depth, kept share, reads and bases kept.
fn write_fastx_record(
    w: &mut BufWriter<fs::File>,
    id: &[u8],
    seq: &[u8],
    qual: Option<&[u8]>,
) -> Result<()> {
    if let Some(q) = qual.filter(|q| q.len() == seq.len()) {
        w.write_all(b"@")?;
        w.write_all(id)?;
        w.write_all(b"\n")?;
        w.write_all(seq)?;
        w.write_all(b"\n+\n")?;
        w.write_all(q)?;
        w.write_all(b"\n")?;
    } else {
        w.write_all(b">")?;
        w.write_all(id)?;
        w.write_all(b"\n")?;
        w.write_all(seq)?;
        w.write_all(b"\n")?;
    }
    Ok(())
}

fn cap_depth(reads: &[PathBuf], out: &Path, salt: u64) -> Result<(Vec<PathBuf>, Value)> {
    let mut records: Vec<(Vec<u8>, Vec<u8>, Option<Vec<u8>>)> = Vec::new();
    for path in reads {
        if !path.is_file() || fs::metadata(path)?.len() == 0 {
            continue;
        }
        crate::recruit::for_each_rec(path, |rec| {
            records.push((rec.id, rec.seq, rec.qual));
        })?;
    }
    let depths = read_depths(&records.iter().map(|r| r.1.as_slice()).collect::<Vec<_>>());
    let depth = depth_quantile(&depths, 0.5);
    let Some(d) = depth.filter(|&d| d > CORRECT_MAX_DEPTH) else {
        return Ok((reads.to_vec(), json!({"depth": depth, "kept_fraction": 1.0})));
    };
    let fraction = CORRECT_MAX_DEPTH / d;
    // `run ont --read-set target` caps before anything else has made the correction folder
    if let Some(dir) = out.parent() {
        fs::create_dir_all(dir)?;
    }
    let mut w = BufWriter::new(fs::File::create(out)?);
    let (mut n, mut b) = (0u64, 0u64);
    for (id, seq, qual) in records.iter().filter(|r| read_share(&r.0, salt) < fraction) {
        write_fastx_record(&mut w, id, seq, qual.as_deref())?;
        n += 1;
        b += seq.len() as u64;
    }
    w.flush()?;
    log(&format!(
        "reads for {} at depth {d:.0}: correcting {n} of {} ({:.1} Mb), down to depth \
         {CORRECT_MAX_DEPTH:.0}",
        out.display(),
        records.len(),
        b as f64 / 1e6
    ));
    Ok((
        vec![out.to_path_buf()],
        json!({"depth": d, "kept_fraction": fraction, "kept_reads": n, "kept_bases": b,
               "reads": path_str(out)}),
    ))
}

/// `find recruit --preset ont`'s classification: a read is the target's when the bait covers
/// 80% of it, partly covered (a NUMT candidate) when it holds a 2 kb block of the bait.
fn ont_params() -> ClassifyParams {
    ClassifyParams {
        min_cover: 0.8,
        min_hit_frac: 0.05,
        numt_min_block: 2000,
        block_gap: 200,
        min_distinct: 0.0,
    }
}

/// The partly covered reads of a recruitment pass (`numt_candidates.fastq`) by the target whose
/// bait they fit best, the bait rebuilt as the pass had it; at most about `max_bases` (0: all)
/// per target. A read that reaches past its target's bait is recruited only when the bait covers
/// 80% of it, so the reads past a bait's edge are a fifth of the depth there and fall below the
/// solid count of correction and of `--seed-reads`; with the partly covered ones the depth
/// stays for kilobases. NUMT reads are among them, at nuclear depth, so their flanks stay out of
/// the solid k-mers. Per target: (FASTA, reads, bases).
fn split_partial(
    dir: &Path,
    k: usize,
    seeds: &BTreeMap<String, Vec<PathBuf>>,
    seed_reads: &[(String, PathBuf)],
    max_bases: u64,
    salt: u64,
) -> Result<BTreeMap<String, (PathBuf, u64, u64)>> {
    let mut bait = BaitSet::new(k);
    for (name, paths) in seeds {
        let t = bait.target(name)?;
        for p in paths {
            bait.add_seed_file(t, p)?;
        }
    }
    for (name, p) in seed_reads {
        let t = bait.target(name)?;
        bait.add_seed_reads(t, p, 3)?;
    }
    bait.prepare_lookup();
    let params = ont_params();
    let path = dir.join("numt_candidates.fastq");
    let mut calls: Vec<(Vec<u8>, Vec<u8>, Option<Vec<u8>>, usize)> = Vec::new();
    let mut bases = vec![0u64; bait.targets.len()];
    let meta = fs::metadata(&path).with_context(|| format!("cannot access {}", path.display()))?;
    anyhow::ensure!(meta.is_file(), "not a regular file: {}", path.display());
    if meta.len() > 0 {
        let mut reader = needletail::parse_fastx_file(&path)
            .with_context(|| format!("cannot read {}", path.display()))?;
        while let Some(record) = reader.next() {
            let record =
                record.with_context(|| format!("malformed record in {}", path.display()))?;
            let seq = record.seq();
            if let (Call::NumtCandidate(t), _) = classify(&seq, &bait, &params) {
                calls.push((
                    record.id().to_vec(),
                    seq.into_owned(),
                    record.qual().map(|q| q.to_vec()),
                    t,
                ));
            }
        }
    }
    // only reads at the depth of what their target recruited (see `at_target_depth`)
    let mut recruited_all: Vec<Vec<Vec<u8>>> = Vec::with_capacity(bait.targets.len());
    for (t, target) in bait.targets.iter().enumerate() {
        let recruited = read_seqs(&dir.join(format!("{}.fastq", target.name)))?;
        let idx: Vec<usize> = (0..calls.len()).filter(|&i| calls[i].3 == t).collect();
        let keep = at_target_depth(
            &recruited,
            &idx.iter().map(|&i| calls[i].1.as_slice()).collect::<Vec<_>>(),
        );
        let dropped = keep.iter().filter(|&&k| !k).count();
        if dropped > 0 {
            log(&format!(
                "{dropped} of {} partly covered {} reads are below a quarter of its recruited \
                 reads' depth (another organelle's or nuclear) and are left out",
                idx.len(),
                target.name
            ));
        }
        for (&i, k) in idx.iter().zip(keep) {
            if !k {
                calls[i].3 = usize::MAX;
            }
        }
        recruited_all.push(recruited);
    }
    if bait.targets.len() > 1 {
        let report_path = dir.join("recruit.json");
        let report_json = report_path.is_file().then(|| read_json(&report_path)).transpose()?;
        let fractions: Vec<f64> = bait
            .targets
            .iter()
            .map(|target| {
                report_json
                    .as_ref()
                    .and_then(|r| target_output(r, &target.name).ok())
                    .and_then(|o| o["sampling_fraction"].as_f64())
                    .unwrap_or(1.0)
            })
            .collect();
        let active: Vec<usize> = (0..calls.len()).filter(|&i| calls[i].3 != usize::MAX).collect();
        let active_refs: Vec<(&[u8], usize)> = active
            .iter()
            .map(|&i| (calls[i].1.as_slice(), calls[i].3))
            .collect();
        let keep_cross = filter_cross_organelle_partials(
            &recruited_all,
            &fractions,
            &active_refs,
            params.numt_min_block,
        );
        let mut dropped_by_target = vec![0usize; bait.targets.len()];
        let mut before_by_target = vec![0usize; bait.targets.len()];
        for (&i, keep) in active.iter().zip(keep_cross) {
            let t = calls[i].3;
            before_by_target[t] += 1;
            if !keep {
                dropped_by_target[t] += 1;
                calls[i].3 = usize::MAX;
            }
        }
        for (t, target) in bait.targets.iter().enumerate() {
            if dropped_by_target[t] > 0 {
                log(&format!(
                    "{} of {} partly covered {} reads stand at another deeper organelle's \
                     normalized depth without a {} bp {}-unique flank and are left out",
                    dropped_by_target[t],
                    before_by_target[t],
                    target.name,
                    params.numt_min_block,
                    target.name
                ));
            }
        }
    }
    for c in &calls {
        if c.3 != usize::MAX {
            bases[c.3] += c.1.len() as u64;
        }
    }
    let mut out = BTreeMap::new();
    for (t, target) in bait.targets.iter().enumerate() {
        let fraction = if max_bases > 0 && bases[t] > max_bases {
            max_bases as f64 / bases[t] as f64
        } else {
            1.0
        };
        let fasta = dir.join(format!("{}.partial.fasta", target.name));
        let mut w = BufWriter::new(fs::File::create(&fasta)?);
        let (mut n, mut b) = (0u64, 0u64);
        for (id, seq, qual, _) in calls.iter().filter(|c| c.3 == t) {
            if fraction < 1.0 && read_share(id, salt) >= fraction {
                continue;
            }
            write_fastx_record(&mut w, id, seq, qual.as_deref())?;
            n += 1;
            b += seq.len() as u64;
        }
        w.flush()?;
        if n > 0 {
            out.insert(target.name.clone(), (fasta, n, b));
        }
    }
    Ok(out)
}

/// The output of `organelle` in a `recruit.json`.
fn target_output<'a>(report: &'a Value, organelle: &str) -> Result<&'a Value> {
    let index = report["targets"]
        .as_array()
        .into_iter()
        .flatten()
        .position(|t| t["name"] == organelle)
        .with_context(|| format!("recruit.json has no target {organelle}"))?;
    Ok(&report["outputs"][index])
}

/// Whether a recruitment round is worth using after one that recruited `previous` bases.
pub fn round_grew(previous: Option<u64>, now: u64) -> bool {
    match previous {
        None => true,
        Some(p) => now > p && now as f64 >= p as f64 * (1.0 + ROUND_GROWTH),
    }
}

fn seed_spec(name: &str, path: &Path) -> OsString {
    let mut spec = OsString::from(format!("{name}="));
    spec.push(path);
    spec
}

/// `assemble correct` of `reads` into `fasta` (and its report beside it): self-correction, or
/// hybrid with `short` (then the self rounds of `assemble correct`'s default); the reads down to
/// `CORRECT_MAX_DEPTH` first (see `cap_depth`).
fn correct(
    reads: &[PathBuf],
    short: Option<&Path>,
    chemistry: crate::r10_correct::OntChemistry,
    fasta: &Path,
    salt: u64,
) -> Result<Value> {
    let json = fasta.with_extension("json");
    let (reads, depth_cap) = cap_depth(reads, &fasta.with_extension("capped.fasta"), salt)?;
    let mut argv = Argv::new("assemble correct")
        .each("--reads", &reads)
        .opt("--out", fasta)
        .opt("--out-json", &json);
    if let Some(s) = short {
        argv = argv.opt("--short-reads", s);
    } else {
        argv = argv.opt("--chemistry", chemistry.as_arg());
    }
    argv.run()?;
    let report = read_json(&json)?;
    let last = report["rounds"].as_array().and_then(|r| r.last()).cloned();
    let count = last.as_ref().and_then(|l| l["reads_out"].as_u64())
        .ok_or_else(|| anyhow::anyhow!("correction report has no valid final reads_out: {}", json.display()))?;
    if short.is_none() { require_self_correction(count, fasta)?; }
    let mut out = json!({
        "reads": path_str(fasta),
        "report": path_str(&json),
        "chemistry": report["chemistry"],
        "hybrid": short.is_some(),
        "depth_cap": depth_cap,
        "rounds": report["rounds"].as_array().map(|r| r.iter().map(|x| json!({
            "k": x["k"],
            "hybrid": x["hybrid"],
            "solid_peak": x["solid_peak"],
            "reads_out": x["reads_out"],
        })).collect::<Vec<_>>()),
        "reads_out": last.as_ref().map(|l| l["reads_out"].clone()),
        "bases_out": last.as_ref().map(|l| l["stats"]["bases_out"].clone()),
        "seconds": report["elapsed_seconds"],
    });
    if !report["r10"].is_null() {
        out["r10"] = report["r10"].clone();
    }
    Ok(out)
}

fn correction_output_usable(count: u64, fasta: &Path) -> Result<bool> {
    let bytes = std::fs::metadata(fasta)?.len();
    Ok(count > 0 && bytes > 0)
}

fn require_self_correction(count: u64, fasta: &Path) -> Result<()> {
    if !correction_output_usable(count, fasta)? {
        anyhow::bail!("self-correction produced no reads: {}; assembly cannot continue", fasta.display());
    }
    Ok(())
}

fn use_hybrid_correction(count: u64, fasta: &Path, self_corrected: &Path, report: &mut Value) -> Result<bool> {
    if correction_output_usable(count, fasta)? { return Ok(true); }
    log("hybrid correction produced no reads; retaining self-corrected reads");
    report["fallback"] = json!({"reason": "hybrid correction produced no reads", "reads": path_str(self_corrected)});
    Ok(false)
}

fn path_str(p: &Path) -> String {
    p.display().to_string()
}

/// Per organelle: its reads so far, or why it has none.
struct Reads {
    /// Uncorrected long reads of the organelle (the last recruitment round that grew).
    raw: Vec<PathBuf>,
    /// Bait of the next round: those reads and the partly covered ones, corrected.
    bait: Option<PathBuf>,
    /// The reads that are assembled: `raw`, corrected (hybrid in `run hybrid`).
    corrected: PathBuf,
    kept: Option<(u64, u64)>,
    /// A used round reached `--recruit-max-bases`: later rounds cannot add to what is kept, so
    /// they are not used for it and it gets no bait of reads.
    capped: bool,
    rounds: Vec<Value>,
    corrections: Vec<Value>,
    short: Value,
}

fn run(long: &[PathBuf], short: Option<Short>, o: &LongOptions) -> Result<Value> {
    let t0 = Instant::now();
    validate(long, short.as_ref(), o)?;
    if let Some(n) = o.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n as usize)
            .build_global()?;
    }
    let pipeline_name = if short.is_some() { "hybrid" } else { "ont" };
    let self_chemistry = if short.is_some() {
        crate::r10_correct::OntChemistry::Legacy
    } else {
        o.chemistry
    };
    let organelles: &[&'static str] = match o.organelle {
        OrganelleChoice::Mitochondrion => &["mitochondrion"],
        OrganelleChoice::Plastid => &["plastid"],
        OrganelleChoice::Both => &["mitochondrion", "plastid"],
    };
    let single = (organelles.len() == 1).then(|| organelles[0]);
    let mut seeds = seed_map(&o.seeds, single)?;
    fs::create_dir_all(&o.out)?;
    // with --read-set target too: the seeds are also the assembly's depth reference
    // (`pipeline::run` passes the first, as `choose_assemble_k` does). Without one the node depth
    // of corrected noisy reads is their error nodes' (Nipponbare ONT plastid: 8x where the
    // sequence stood at 430x, min count 1, 833 unitigs, 46% explained)
    let builtin = crate::seeddb::fill_missing(
        &mut seeds,
        &o.out.join("seeds"),
        single,
        o.single_seed,
        o.lineage,
    )?;
    {
        if !builtin.is_empty() {
            log(&format!(
                "no --seeds for {}: using the built-in {} seed library {}",
                builtin.iter().cloned().collect::<Vec<_>>().join(" and "),
                o.lineage.name(),
                o.lineage.version()
            ));
        }
    }

    let mut reads: BTreeMap<&'static str, Result<Reads>> = BTreeMap::new();
    let mut files = BTreeMap::new();
    match o.read_set {
        ReadSet::Target => {
            let organelle = organelles[0];
            log(&format!("correcting the supplied {organelle} long reads"));
            let fasta = o.out.join(format!("correct/{organelle}.fasta"));
            let mut c = correct(long, None, self_chemistry, &fasta, o.salt)?;
            c["purpose"] = "assembly".into();
            reads.insert(
                organelle,
                Ok(Reads {
                    raw: long.to_vec(),
                    bait: None,
                    corrected: fasta,
                    kept: None,
                    capped: false,
                    rounds: Vec::new(),
                    corrections: vec![c],
                    short: Value::Null,
                }),
            );
        }
        ReadSet::WholeGenome => {
            for round in 1..=o.recruit_rounds {
                let dir = o.out.join("recruit").join(format!("round{round}"));
                // the corrected reads that grew in the round before, as bait
                let seed_reads: Vec<(String, PathBuf)> = reads
                    .iter()
                    .filter_map(|(name, r)| {
                        r.as_ref().ok()?.bait.clone().map(|b| (name.to_string(), b))
                    })
                    .collect();
                log(&format!(
                    "long-read recruitment round {round}: {}",
                    if seed_reads.is_empty() {
                        "seeds"
                    } else {
                        "seeds and the corrected reads of the round before"
                    }
                ));
                let mut argv = Argv::new("find recruit")
                    .opt("--out", &dir)
                    .opt("--preset", "ont")
                    .opt("--k", o.recruit_k.to_string())
                    .opt("--iterations", "1")
                    .opt("--salt", o.salt.to_string())
                    .each("--reads", long);
                if o.recruit_max_bases > 0 {
                    argv = argv.opt("--max-bases", o.recruit_max_bases.to_string());
                }
                if let Some(m) = o.max_read_methylation {
                    argv = argv.opt("--max-read-methylation", m.to_string());
                }
                for (name, paths) in &seeds {
                    for p in paths {
                        argv = argv.opt("--seed", seed_spec(name, p));
                    }
                }
                for (name, p) in &seed_reads {
                    argv = argv.opt("--seed-reads", seed_spec(name, p));
                }
                argv.run()?;
                files.insert(format!("recruitment_round{round}"), path_str(&dir.join("recruit.json")));
                let report = read_json(&dir.join("recruit.json"))?;
                let mut grown = Vec::new();
                let mut capped = Vec::new();
                for &organelle in organelles {
                    if matches!(reads.get(organelle), Some(Err(_))) {
                        continue;
                    }
                    let out = target_output(&report, organelle)?;
                    let n = |key: &str| out[key].as_u64().unwrap_or(0);
                    let before = reads.get(organelle).and_then(|r| r.as_ref().ok());
                    let previous = before
                        .and_then(|r| r.rounds.iter().rev().find(|x| x["used"] == true))
                        .and_then(|r| r["recruited_bases"].as_u64());
                    let was_capped = before.is_some_and(|r| r.capped);
                    let grew = !was_capped && round_grew(previous, n("recruited_bases"));
                    let at_cap = grew && out["sampled_by"] == "max_bases";
                    let entry = json!({
                        "round": round,
                        "bait": if seed_reads.iter().any(|(s, _)| s == organelle) {
                            "seeds + corrected reads"
                        } else {
                            "seeds"
                        },
                        "recruited_reads": n("recruited_reads"),
                        "recruited_bases": n("recruited_bases"),
                        "kept_reads": n("kept_reads"),
                        "kept_bases": n("kept_bases"),
                        "sampled_by": out["sampled_by"],
                        "used": grew,
                        "capped_before": was_capped,
                        "seconds": report["elapsed_seconds"],
                    });
                    let raw = vec![dir.join(format!("{organelle}.fastq"))];
                    let kept = Some((n("kept_reads"), n("kept_bases")));
                    match reads.get_mut(organelle) {
                        Some(Ok(r)) => {
                            r.rounds.push(entry);
                            if grew {
                                (r.raw, r.kept, r.capped) = (raw, kept, at_cap);
                            }
                        }
                        // kept even with no read covered by 80%: with a multi-species seed library
                        // a mitochondrial ONT read rarely is (ZM4: none in 12 Gb), and the partly
                        // covered reads alone give the next round's bait
                        _ => {
                            reads.insert(
                                organelle,
                                Ok(Reads {
                                    raw,
                                    bait: None,
                                    corrected: PathBuf::new(),
                                    kept,
                                    capped: at_cap,
                                    rounds: vec![entry],
                                    corrections: Vec::new(),
                                    short: Value::Null,
                                }),
                            );
                        }
                    }
                    if at_cap {
                        capped.push(organelle);
                    } else if grew {
                        grown.push(organelle);
                    }
                }
                if !capped.is_empty() && round < o.recruit_rounds {
                    log(&format!(
                        "round {round} reached --recruit-max-bases for {}: its reads are kept and \
                         later rounds are not used for it",
                        capped.join(" and ")
                    ));
                }
                if grown.is_empty() {
                    if round < o.recruit_rounds {
                        log(&format!(
                            "round {round}: no organelle below --recruit-max-bases gained {:.0}% \
                             of its bases; the reads of each one's last used round are assembled",
                            ROUND_GROWTH * 100.0
                        ));
                    }
                    break;
                }
                if round == o.recruit_rounds {
                    break;
                }
                // bait of the next round: this round's reads and the partly covered ones, which
                // reach past the bait at organelle depth, corrected together
                let partial = split_partial(
                    &dir,
                    o.recruit_k,
                    &seeds,
                    &seed_reads,
                    o.recruit_max_bases,
                    o.salt,
                )?;
                for organelle in grown {
                    let Some(Ok(r)) = reads.get_mut(organelle) else { continue };
                    let mut inputs =
                        if r.kept.is_some_and(|k| k.1 > 0) { r.raw.clone() } else { Vec::new() };
                    if let Some((p, n, _)) = partial.get(organelle) {
                        log(&format!(
                            "self-correcting the {organelle} long reads of round {round} with \
                             {n} partly covered reads, as bait for the next round"
                        ));
                        inputs.push(p.clone());
                    }
                    if inputs.is_empty() {
                        continue;
                    }
                    let fasta = o.out.join(format!("correct/round{round}/{organelle}.bait.fasta"));
                    let result = (|| -> Result<(Value, PathBuf)> {
                        let mut c = correct(&inputs, None, self_chemistry, &fasta, o.salt)?;
                        c["purpose"] = "bait".into();
                        c["partial_reads"] = partial.get(organelle).map(|x| x.1).into();
                        let masked = fasta.with_extension("masked.fasta");
                        c["tandem_masked_bases"] = mask_tandem_fasta(&fasta, &masked)?.into();
                        Ok((c, masked))
                    })();
                    match result {
                        Ok((c, masked)) => { r.bait = Some(masked); r.corrections.push(c); }
                        Err(e) => { reads.insert(organelle, Err(e.context("bait correction"))); }
                    }
                }
                // A capped organelle (e.g. plastid reaching --recruit-max-bases in round 1) keeps
                // this round's uniform sample for assembly, but must also self-correct its own
                // sample-specific bait before later rounds run: otherwise round 2 competes a
                // sample-specific mitochondrial bait (holding exact sample MTPTs) against only
                // multi-species plastid seeds, causing thousands of pure plastid reads over MTPTs
                // to defect into mitochondrion and overwrite mitochondrial MTPT alleles.
                for &organelle in &capped {
                    let Some(Ok(r)) = reads.get_mut(organelle) else { continue };
                    if r.bait.is_some() || r.kept.map_or(0, |k| k.1) < MIN_RECRUITED_BASES {
                        continue;
                    }
                    log(&format!(
                        "self-correcting the {organelle} long reads to assemble, and as competitive bait for the next round"
                    ));
                    let fasta = o.out.join(format!("correct/{organelle}.fasta"));
                    let result = (|| -> Result<(Value, PathBuf)> {
                        let mut c = correct(&r.raw, None, self_chemistry, &fasta, o.salt)?;
                        c["purpose"] = "assembly".into();
                        let masked = fasta.with_extension("masked.fasta");
                        c["tandem_masked_bases"] = mask_tandem_fasta(&fasta, &masked)?.into();
                        Ok((c, masked))
                    })();
                    match result {
                        Ok((c, masked)) => {
                            r.corrected = fasta;
                            r.bait = Some(masked);
                            r.corrections.push(c);
                        }
                        Err(e) => {
                            reads.insert(organelle, Err(e.context("self correction for assembly")));
                        }
                    }
                }
            }
            for (&organelle, entry) in reads.iter_mut() {
                let Ok(r) = entry else { continue };
                let bases = r.kept.map_or(0, |k| k.1);
                if bases < MIN_RECRUITED_BASES {
                    *entry = Err(anyhow::anyhow!(
                        "only {bases} bases of {organelle} long reads were recruited (at least \
                         {MIN_RECRUITED_BASES} are needed) after {} round(s): the seeds may be \
                         too distant for noisy reads; give closer seeds or more --recruit-rounds",
                        r.rounds.len()
                    ));
                    continue;
                }
                if !r.corrected.as_os_str().is_empty() {
                    continue;
                }
                log(&format!("self-correcting the {organelle} long reads to assemble"));
                let fasta = o.out.join(format!("correct/{organelle}.fasta"));
                match correct(&r.raw, None, self_chemistry, &fasta, o.salt) {
                    Ok(mut c) => { c["purpose"] = "assembly".into(); r.corrected = fasta; r.corrections.push(c); }
                    Err(e) => { *entry = Err(e.context("self correction for assembly")); }
                }
            }
        }
    }

    if let Some(s) = &short {
        let dir = o.out.join("recruit/short");
        log("recruiting short reads with the corrected long reads");
        let mut argv = Argv::new("find recruit")
            .opt("--out", &dir)
            .opt("--preset", "sr")
            .opt("--iterations", "1")
            .arg("--mates-together")
            .opt("--salt", o.salt.to_string())
            .each("--reads", s.reads);
        if s.max_bases > 0 {
            argv = argv.opt("--max-bases", s.max_bases.to_string());
        }
        let mut baited = 0;
        for (name, entry) in &mut reads {
            if let Ok(r) = entry {
                let masked = r.corrected.with_extension("masked.fasta");
                match mask_tandem_fasta(&r.corrected, &masked) {
                    Ok(_) => { argv = argv.opt("--seed-reads", seed_spec(name, &masked)); baited += 1; }
                    Err(e) => { *entry = Err(e.context("masking corrected reads for short-read bait")); }
                }
            }
        }
        // an organelle without corrected reads (the competing one) claims its short reads by its seeds
        for (name, paths) in &seeds {
            if !matches!(reads.get(name.as_str()), Some(Ok(_))) {
                for p in paths {
                    argv = argv.opt("--seed", seed_spec(name, p));
                }
            }
        }
        if baited > 0 {
            argv.run()?;
            files.insert("short_read_recruitment".into(), path_str(&dir.join("recruit.json")));
            let report = read_json(&dir.join("recruit.json"))?;
            for &organelle in organelles {
                let Some(Ok(r)) = reads.get_mut(organelle) else { continue };
                let out = target_output(&report, organelle)?;
                let kept = out["kept_bases"].as_u64().unwrap_or(0);
                r.short = json!({
                    "reads": path_str(&dir.join(format!("{organelle}.fastq"))),
                    "recruited_reads": out["recruited_reads"],
                    "recruited_bases": out["recruited_bases"],
                    "kept_reads": out["kept_reads"],
                    "kept_bases": kept,
                    "sampled_by": out["sampled_by"],
                    "seconds": report["elapsed_seconds"],
                });
                if kept < MIN_RECRUITED_BASES {
                    log(&format!(
                        "only {kept} short-read bases of the {organelle}: its long reads keep their self-correction"
                    ));
                    continue;
                }
                log(&format!("correcting the {organelle} long reads with its short reads"));
                let fasta = o.out.join("correct/hybrid").join(format!("{organelle}.fasta"));
                let result = (|| -> Result<(Value, bool)> {
                    let mut c = correct(
                        &r.raw,
                        Some(&dir.join(format!("{organelle}.fastq"))),
                        crate::r10_correct::OntChemistry::Legacy,
                        &fasta,
                        o.salt,
                    )?;
                    c["purpose"] = "assembly".into();
                    let n = c["reads_out"].as_u64().ok_or_else(|| anyhow::anyhow!("invalid hybrid reads_out"))?;
                    let use_corrected = use_hybrid_correction(n, &fasta, &r.corrected, &mut c)?;
                    Ok((c, use_corrected))
                })();
                match result {
                    Ok((c, use_corrected)) => { if use_corrected { r.corrected = fasta; } r.corrections.push(c); }
                    Err(e) => { reads.insert(organelle, Err(e.context("hybrid correction"))); }
                }
            }
        }
    }

    let mut results = Map::new();
    let mut failed = Vec::new();
    for &organelle in organelles {
        let dir = if single.is_some() { o.out.clone() } else { o.out.join(organelle) };
        let result = reads
            .remove(organelle)
            .unwrap_or_else(|| Err(anyhow::anyhow!("no reads of the {organelle}")))
            .and_then(|r| {
                assemble(o, long, organelle, &r, &dir, &seeds, &builtin, files.clone(), pipeline_name, t0)
            });
        match result {
            Ok(report) => {
                results.insert(organelle.into(), report);
            }
            Err(e) if single.is_some() => return Err(e),
            Err(e) => {
                let message = format!("{e:#}");
                log(&format!("{organelle} failed: {message}"));
                results.insert(organelle.into(), json!({"error": message}));
                failed.push(organelle);
            }
        }
    }
    if let Some(organelle) = single {
        return Ok(results.remove(organelle).expect("the one organelle"));
    }
    let report = json!({
        "backend": "ovasm",
        "pipeline": pipeline_name,
        "organelles": results,
        "cross_organelle_candidates": crate::retained::cross_organelle(&o.out, &o.out, crate::kmer::MAX_K)?,
        "failed": failed,
        "elapsed_seconds": t0.elapsed().as_secs_f64(),
        "peak_rss_mb": peak_rss_mb(),
    });
    write_json(&o.out.join("summary.json"), &report)?;
    if !failed.is_empty() {
        bail!(
            "{} failed (see {}); the other organelle's results are in its directory",
            failed.join(" and "),
            o.out.join("summary.json").display()
        );
    }
    Ok(report)
}

/// k of the graph of corrected `reads`: `ASSEMBLE_K_HIGH`, or `ASSEMBLE_K_LOW` when its graph holds
/// `ASSEMBLE_K_GROWTH` times the sequence of the higher k's (see there). Each probe is one
/// `assemble graph` pass over the reads, written beside the graph and removed.
fn choose_assemble_k(
    reads: &Path,
    dir: &Path,
    depth_reference: Option<&Path>,
    is_r10: bool,
) -> Result<usize> {
    let length = |k: usize| -> Result<f64> {
        let (gfa, json) = (
            dir.join(format!("probe-k{k}.gfa")),
            dir.join(format!("probe-k{k}.json")),
        );
        let mut probe = Argv::new("assemble graph")
            .opt("--reads", reads)
            .opt("--k", k.to_string())
            .opt("--out", &gfa)
            .opt("--out-json", &json);
        // the depth the assembly itself uses: without it the nodes' plain median counts the
        // error nodes of noisy reads (21x on Nipponbare ONT against 370x on the seeds' own
        // nodes), the minimum count falls to 3 and both graphs fill with noise
        if let Some(reference) = depth_reference {
            probe = probe.opt("--depth-reference", reference);
        }
        probe.run()?;
        let total = read_json(&json)?["total_length"].as_f64().unwrap_or(0.0);
        for f in [&gfa, &json, &gfa.with_extension("cleaning.fa")] {
            let _ = std::fs::remove_file(f);
        }
        Ok(total)
    };
    let (high, low) = (length(ASSEMBLE_K_HIGH)?, length(ASSEMBLE_K_LOW)?);
    let mut k = pick_assemble_k(high, low);
    if is_r10 && k == ASSEMBLE_K_HIGH {
        let r10_len = length(ASSEMBLE_K_R10)?;
        if r10_len > 0.0 && high < ASSEMBLE_K_GROWTH * r10_len {
            k = ASSEMBLE_K_R10;
        }
        log(&format!(
            "graph of {} bp at k {ASSEMBLE_K_R10}, {} bp at k {ASSEMBLE_K_HIGH}, {} bp at k {ASSEMBLE_K_LOW}: starting at k {k}",
            r10_len as u64, high as u64, low as u64
        ));
    } else {
        log(&format!(
            "graph of {} bp at k {ASSEMBLE_K_HIGH}, {} bp at k {ASSEMBLE_K_LOW}: starting at k {k}",
            high as u64, low as u64
        ));
    }
    Ok(k)
}

/// `ASSEMBLE_K_LOW` when its graph (`low` bases) holds `ASSEMBLE_K_GROWTH` times the sequence of
/// the higher k's (`high`), else `ASSEMBLE_K_HIGH`; a higher k with no graph is not preferred.
fn pick_assemble_k(high: f64, low: f64) -> usize {
    if low > 0.0 && (high <= 0.0 || low >= ASSEMBLE_K_GROWTH * high) {
        ASSEMBLE_K_LOW
    } else {
        ASSEMBLE_K_HIGH
    }
}

/// The HiFi stages on one organelle's corrected reads, then its summary with the long-read
/// record (and the short-read evidence of `run hybrid`).
#[allow(clippy::too_many_arguments)]
fn assemble(
    o: &LongOptions,
    long: &[PathBuf],
    organelle: &str,
    r: &Reads,
    dir: &Path,
    seeds: &BTreeMap<String, Vec<PathBuf>>,
    builtin: &BTreeSet<String>,
    mut files: BTreeMap<String, String>,
    pipeline_name: &str,
    t0: Instant,
) -> Result<Value> {
    let organelle_choice = match organelle {
        "mitochondrion" => OrganelleChoice::Mitochondrion,
        _ => OrganelleChoice::Plastid,
    };
    let own_seeds: BTreeMap<String, Vec<PathBuf>> = seeds
        .get(organelle)
        .map(|s| (organelle.to_string(), s.clone()))
        .into_iter()
        .collect();
    let used_r10 = pipeline_name == "ont"
        && r.corrections
            .iter()
            .rev()
            .find(|c| c["purpose"] == "assembly")
            .is_some_and(|c| c["chemistry"] == "r10");
    let assemble_k = match o.assemble_k {
        Some(k) => k,
        None => choose_assemble_k(
            &r.corrected,
            dir,
            own_seeds.get(organelle).and_then(|v| v.first()).map(PathBuf::as_path),
            used_r10,
        )?,
    };
    log(&format!("assembling the corrected {organelle} reads from k {assemble_k}"));
    let mut summary = pipeline::run(&pipeline::RunConfig {
        reads: vec![r.corrected.clone()],
        organelle: organelle_choice,
        out: dir.to_path_buf(),
        read_type: ReadType::Hifi,
        read_set: ReadSet::Target,
        seeds: own_seeds,
        builtin_seeds: builtin.clone(),
        seed_source_audits: Default::default(),
        read_admission_manifest: None,
        read_source_exclusions: None,
        seed_lineage: o.lineage,
        recruit: pipeline::RecruitBy::Seeds,
        recruit_depth: None,
        recruit_cache: None,
        auto_recruit_cache: false,
        gated_recruit: false,
        recruit_gate_fold: None,
        recruit_nuclear_depth: None,
        recruit_genome_size: None,
        assemble_k: Some(assemble_k),
        recruit_max_depth: None,
        recruit_max_bases: None,
        recruit_iterations: None,
        read_context: None,
        single_seed: true,
        mode: o.mode,
        replicates: 3,
        max_replicates: crate::replicate::MAX_REPLICATES,
        replicate_depth: None,
        salt: o.salt,
        sample: o.sample.clone(),
        component_panel: None,
        prune_gene_free: false,
        other_organelle_filter: false,
        shallow_seeds: false,
        depth_reference: Default::default(),
        methylation: o.methylation,
        methylation_bam: if !o.methylation_bam.is_empty() {
            o.methylation_bam.clone()
        } else {
            long.iter().filter(|p| crate::methylation::is_bam_path(p)).cloned().collect()
        },
        max_read_methylation: o.max_read_methylation,
    })?;
    // homopolymer runs from the reads (see `hp_polish`): the molecules the stages wrote are kept
    // beside as molecules.unpolished.fasta. With short reads of the sample (`run hybrid`) those
    // judge the runs: they read homopolymers right, and the hybrid correction has already set
    // the runs they cover; the uncorrected long reads would set them back to their own systematic
    // lengths (code audit 11-15: 11 homopolymer differences written into the plastid). A run the
    // short reads do not reach `MIN_OBSERVATIONS` times stays as the hybrid correction left it.
    let molecules = dir.join("molecules.fasta");
    let short_for_polish = r.short["reads"]
        .as_str()
        .filter(|_| r.short["kept_bases"].as_u64().unwrap_or(0) >= MIN_RECRUITED_BASES)
        .map(PathBuf::from);
    let (polish_reads, polish_with) = match short_for_polish {
        Some(sr) => (vec![sr], "short reads"),
        None => (r.raw.clone(), "uncorrected reads"),
    };
    let polish = if molecules.exists() && used_r10 {
        json!({
            "reads_used": "r10 phased consensus",
            "skipped": true,
        })
    } else if molecules.exists() {
        log(&format!("polishing the homopolymer runs of the {organelle} molecules with the {polish_with}"));
        let unpolished = dir.join("molecules.unpolished.fasta");
        fs::rename(&molecules, &unpolished)?;
        let rep = crate::hp_polish::polish(&unpolished, &polish_reads, &molecules, crate::hp_polish::MIN_OBSERVATIONS)?;
        let (added, removed) = rep.molecules.iter().fold((0, 0), |a, m| (a.0 + m.bases_added, a.1 + m.bases_removed));
        log(&format!(
            "{} runs judged from {} reads; {} bases added, {} removed",
            rep.molecules.iter().map(|m| m.runs_judged).sum::<usize>(),
            rep.reads_placed,
            added,
            removed
        ));
        files.insert("unpolished_molecules".into(), path_str(&unpolished));
        let mut v = serde_json::to_value(&rep)?;
        v["reads_used"] = polish_with.into();
        v
    } else {
        Value::Null
    };
    let mut short_evidence = Value::Null;
    if let Some(sr) = r.short["reads"].as_str() {
        let sr = PathBuf::from(sr);
        if r.short["kept_bases"].as_u64().unwrap_or(0) >= MIN_RECRUITED_BASES {
            log(&format!("checking the {organelle} graph against its short reads"));
            let out = dir.join("evidence.short.json");
            evidence_argv(&dir.join("assembly.gfa"), &[sr], &out, ReadType::Sr).run()?;
            let e = read_json(&out)?;
            files.insert("short_read_evidence".into(), path_str(&out));
            short_evidence = json!({
                "links": e["links"].as_array().map_or(0, Vec::len),
                "links_supported": e["links_supported"],
                "adjacency_supported_links": e["adjacency_supported_links"],
            });
        }
    }
    // --read-set target loads them too, for the depth reference only
    let seed_info = json!({
        "source": if builtin.contains(organelle) { "builtin" } else { "custom" },
        "files": seeds.get(organelle).map(|v| v.iter().map(|p| path_str(p)).collect::<Vec<_>>()),
        "used_for": if o.read_set == ReadSet::Target { "depth_reference" } else { "recruitment_and_depth_reference" },
    });
    let s = summary.as_object_mut().context("summary.json is not an object")?;
    s.insert(
        "read_set".into(),
        match o.read_set {
            ReadSet::WholeGenome => "whole_genome",
            ReadSet::Target => "target_reads",
        }
        .into(),
    );
    s.insert("read_type".into(), if pipeline_name == "hybrid" { "ont_hybrid" } else { "ont" }.into());
    s.insert("seed_database".into(), seed_info);
    s.insert("recruited_reads".into(), r.kept.map(|k| k.0).into());
    s.insert("recruited_bases".into(), r.kept.map(|k| k.1).into());
    s.insert(
        "long_read".into(),
        json!({
            "pipeline": pipeline_name,
            "recruit_k": o.recruit_k,
            "rounds": r.rounds,
            "raw_reads": r.raw.iter().map(|p| path_str(p)).collect::<Vec<_>>(),
            "corrections": r.corrections,
            "assembled_reads": path_str(&r.corrected),
            "short_reads": r.short,
        }),
    );
    s.insert("short_read_evidence".into(), short_evidence);
    s.insert("homopolymer_polish".into(), polish);
    s.insert("total_elapsed_seconds".into(), t0.elapsed().as_secs_f64().into());
    s.insert("peak_rss_mb".into(), peak_rss_mb().into());
    if let Some(Value::Object(f)) = s.get_mut("files") {
        for (k, v) in files {
            f.entry(k).or_insert(v.into());
        }
        f.insert("corrected_reads".into(), path_str(&r.corrected).into());
    }
    write_json(&dir.join("summary.json"), &summary)?;
    Ok(summary)
}

#[cfg(test)]
mod tests {
    #[test]
    fn empty_correction_is_not_usable() {
        let p = std::env::temp_dir().join(format!("ovasm-empty-correction-{}", std::process::id()));
        std::fs::write(&p, b"").unwrap();
        assert!(!super::correction_output_usable(0, &p).unwrap());
        assert!(!super::correction_output_usable(1, &p).unwrap());
        std::fs::write(&p, b">r\nACGT\n").unwrap();
        assert!(!super::correction_output_usable(0, &p).unwrap());
        assert!(super::correction_output_usable(1, &p).unwrap());
        assert!(super::require_self_correction(0, &p).unwrap_err().to_string().contains("self-correction produced no reads"));
        let mut report = serde_json::json!({});
        assert!(!super::use_hybrid_correction(0, &p, &p, &mut report).unwrap());
        assert_eq!(report["fallback"]["reason"], "hybrid correction produced no reads");
        assert!(super::use_hybrid_correction(1, &p, &p, &mut serde_json::json!({})).unwrap());
        std::fs::remove_file(&p).unwrap();
        assert!(super::require_self_correction(0, &p).is_err());
        assert!(super::require_self_correction(1, &p).is_err());
        assert!(super::use_hybrid_correction(0, &p, &p, &mut serde_json::json!({})).is_err());
        assert!(super::use_hybrid_correction(1, &p, &p, &mut serde_json::json!({})).is_err());
    }

    #[test]
    fn a_read_share_follows_its_fraction_for_nanopore_and_other_read_names() {
        // names as the SRA writes Nanopore reads, a plain counter, and Illumina-like names
        let formats: [fn(u64) -> String; 3] = [
            |i| {
                let x = i.wrapping_mul(0x9E37_79B9_7F4A_7C15);
                format!(
                    "SRR25241091.{}/1 {:08x}-{:04x}-{:04x}-{:04x}-{:012x} length={}",
                    96 + 7 * i,
                    x >> 32,
                    (x >> 16) & 0xffff,
                    x & 0xffff,
                    (x >> 8) & 0xffff,
                    x >> 4,
                    500 + (x >> 40) % 30_000
                )
            },
            |i| format!("read{i}"),
            |i| format!("A00123:45:HXXXXDSXX:1:{}:{}:{} 1:N:0:ACGT", 1101 + i % 4, i % 30_000, i / 4),
        ];
        for f in formats {
            let ids: Vec<Vec<u8>> = (0..40_000u64).map(|i| f(i).into_bytes()).collect();
            for fraction in [0.12, 0.16, 0.5] {
                let kept = ids.iter().filter(|id| read_share(id, 0) < fraction).count() as f64;
                let share = kept / ids.len() as f64;
                assert!((share - fraction).abs() < 0.02, "{share} kept for {fraction}");
            }
        }
        // the three real names the shortfall was seen on
        for (id, salt) in [
            ("SRR25241091.96/1 38fdbf00-5340-417f-a18c-139e14d2055a length=37207", 0),
            ("SRR25241091.291/1 bf1054c5-a549-4c41-9b8e-7d95d46366d0 length=33228", 0),
        ] {
            let v = read_share(id.as_bytes(), salt);
            assert!((0.0..1.0).contains(&v));
        }
        assert_ne!(read_share(b"a", 0), read_share(b"a", 1), "the salt changes the draw");
    }

    #[test]
    fn the_lower_k_is_chosen_only_when_the_higher_one_lost_sequence() {
        // Nipponbare ONT: 321,855 bp at 201, 404,510 at 101; ZM4 ONT: 299,707 and 300,178
        assert_eq!(pick_assemble_k(321_855.0, 404_510.0), ASSEMBLE_K_LOW);
        assert_eq!(pick_assemble_k(299_707.0, 300_178.0), ASSEMBLE_K_HIGH);
        assert_eq!(pick_assemble_k(100.0, 111.0), ASSEMBLE_K_LOW, "over a tenth more counts");
        assert_eq!(pick_assemble_k(100.0, 109.0), ASSEMBLE_K_HIGH);
        assert_eq!(pick_assemble_k(0.0, 500.0), ASSEMBLE_K_LOW, "recover the nonempty graph");
        assert_eq!(pick_assemble_k(0.0, 0.0), ASSEMBLE_K_HIGH);
        assert_eq!(pick_assemble_k(500.0, 0.0), ASSEMBLE_K_HIGH);
    }

    use super::*;
    use clap::Parser;

    fn tmpdir(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("ovasm-longread-{name}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        d
    }

    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0 >> 11
        }
        fn unit(&mut self) -> f64 {
            self.next() as f64 / (1u64 << 53) as f64
        }
        fn base(&mut self) -> u8 {
            // top bits: the low bits of an LCG repeat with a short period
            b"ACGT"[(self.next() >> 51) as usize]
        }
    }

    fn seq(n: usize, seed: u64) -> Vec<u8> {
        let mut r = Rng(seed);
        (0..n).map(|_| r.base()).collect()
    }

    /// `truth` with substitutions, insertions and deletions at `rate`.
    fn noisy(truth: &[u8], rate: f64, rng: &mut Rng) -> Vec<u8> {
        let mut out = Vec::with_capacity(truth.len() + truth.len() / 10);
        for &b in truth {
            let x = rng.unit();
            if x < rate / 3.0 {
                out.push(rng.base());
            } else if x < 2.0 * rate / 3.0 {
                out.push(b);
                out.push(rng.base());
            } else if x >= rate {
                out.push(b);
            }
        }
        out
    }

    /// Reads of `len` from the circular `genome` every `step` bases, both strands, with errors.
    fn long_reads(genome: &[u8], len: usize, step: usize, rate: f64, tag: &str, rng: &mut Rng) -> Vec<(String, Vec<u8>)> {
        let circ = [genome, &genome[..len]].concat();
        (0..genome.len())
            .step_by(step)
            .enumerate()
            .map(|(i, st)| {
                let r = noisy(&circ[st..st + len], rate, rng);
                let r = if i % 2 == 1 { crate::evidence::revcomp(&r) } else { r };
                (format!("{tag}{i}"), r)
            })
            .collect()
    }

    fn fasta(path: &Path, records: &[(String, Vec<u8>)]) {
        let text: String = records
            .iter()
            .map(|(id, s)| format!(">{id}\n{}\n", std::str::from_utf8(s).unwrap()))
            .collect();
        fs::write(path, text).unwrap();
    }

    /// A 30 kb mitochondrion at 60x, a 20 kb plastid at 150x and 400 kb of nucleus at 1x, as
    /// 8 kb reads with `rate` errors; the mitochondrial seed lacks 10,000-14,000. Round one
    /// recruits the reads the seed covers by 80%, reaching 1.6 kb into the gap from each side,
    /// and leaves 800 bp without a read; round two's bait reaches that far, and the reads across
    /// what is left (covered 86%) come in.
    struct Sample {
        mito: Vec<u8>,
        plastid: Vec<u8>,
        long: PathBuf,
        mt_seed: PathBuf,
        pt_seed: PathBuf,
    }

    fn sample(dir: &Path, rate: f64) -> Sample {
        let (mito, plastid, nuclear) = (seq(30_000, 11), seq(20_000, 22), seq(400_000, 33));
        let mut rng = Rng(44);
        let mut reads = long_reads(&mito, 8000, 130, rate, "m", &mut rng);
        reads.extend(long_reads(&plastid, 8000, 53, rate, "p", &mut rng));
        reads.extend((0..nuclear.len() - 8000).step_by(8000).enumerate().map(|(i, st)| {
            (format!("n{i}"), noisy(&nuclear[st..st + 8000], rate, &mut rng))
        }));
        let long = dir.join("long.fa");
        fasta(&long, &reads);
        let (mt_seed, pt_seed) = (dir.join("mt.fa"), dir.join("pt.fa"));
        fasta(
            &mt_seed,
            &[
                ("mt_a".into(), mito[..10_000].to_vec()),
                ("mt_b".into(), mito[14_000..].to_vec()),
            ],
        );
        fasta(&pt_seed, &[("pt".into(), plastid.clone())]);
        Sample { mito, plastid, long, mt_seed, pt_seed }
    }

    fn options(args: &[&str]) -> LongOptions {
        #[derive(Parser)]
        struct Wrap {
            #[command(flatten)]
            o: LongOptions,
        }
        Wrap::try_parse_from(std::iter::once("t").chain(args.iter().copied())).unwrap().o
    }

    fn molecules(dir: &Path) -> Vec<Vec<u8>> {
        let text = fs::read_to_string(dir.join("molecules.fasta")).unwrap();
        text.split('>')
            .filter(|r| !r.is_empty())
            .map(|r| r.lines().skip(1).collect::<String>().into_bytes())
            .collect()
    }

    /// Share of the 31-mers of `m` that are in the circular `genome` (either strand).
    fn contained(m: &[u8], genome: &[u8]) -> f64 {
        let circ = [genome, &genome[..31]].concat();
        let mut set = rustc_hash::FxHashSet::default();
        for s in [circ.clone(), crate::evidence::revcomp(&circ)] {
            for w in s.windows(31) {
                set.insert(w.to_vec());
            }
        }
        let n = m.windows(31).count();
        m.windows(31).filter(|w| set.contains(*w)).count() as f64 / n.max(1) as f64
    }

    fn assert_genome(dir: &Path, genome: &[u8], what: &str) {
        let s = read_json(&dir.join("summary.json")).unwrap();
        assert_eq!(s["molecules"], json!(1), "{what}: {s}");
        assert_eq!(s["circular"], json!([true]), "{what}");
        let len = s["lengths"][0].as_u64().unwrap() as i64;
        assert!((len - genome.len() as i64).abs() <= 20, "{what}: {len} bp");
        let m = &molecules(dir)[0];
        assert!(contained(m, genome) > 0.99, "{what}: {:.4}", contained(m, genome));
    }

    #[test]
    fn partly_covered_reads_at_another_depth_are_left_out() {
        // recruited: a 20 kb genome at 50x; partial: 10 reads of it, 10 of a genome seen once
        let (a, other) = (seq(20_000, 5), seq(60_000, 6));
        let circ = [a.as_slice(), &a[..4000]].concat();
        let recruited: Vec<Vec<u8>> = (0..20_000).step_by(80).map(|s| circ[s..s + 4000].to_vec()).collect();
        let mut partial: Vec<&[u8]> = (0..10).map(|i| &circ[i * 1500..i * 1500 + 4000]).collect();
        partial.extend((0..10).map(|i| &other[i * 6000..i * 6000 + 4000]));
        let keep = at_target_depth(&recruited, &partial);
        assert_eq!(keep, [vec![true; 10], vec![false; 10]].concat());
    }

    #[test]
    fn with_nothing_recruited_the_partly_covered_reads_are_compared_with_one_another() {
        // 40 reads of a 20 kb genome every 500 bp (8x among themselves), and 10 NUMT reads:
        // 1.5 kb of it in 2.5 kb of a single-copy nuclear flank
        let (a, other) = (seq(20_000, 5), seq(60_000, 6));
        let circ = [a.as_slice(), &a[..4000]].concat();
        let mut partial: Vec<Vec<u8>> = (0..40).map(|i| circ[i * 500..i * 500 + 4000].to_vec()).collect();
        partial.extend((0..10).map(|i| {
            [&a[i * 1500..i * 1500 + 1500], &other[i * 6000..i * 6000 + 2500]].concat()
        }));
        let refs: Vec<&[u8]> = partial.iter().map(Vec::as_slice).collect();
        assert_eq!(at_target_depth(&[], &refs), [vec![true; 40], vec![false; 10]].concat());
        // all of one depth: all kept
        assert_eq!(at_target_depth(&[], &refs[..40]), vec![true; 40]);
        assert_eq!(at_target_depth(&[], &[]), Vec::<bool>::new());
    }

    #[test]
    fn short_tandem_repeats_are_masked_and_the_rest_is_kept() {
        let flank = seq(200, 9);
        let ga = b"GA".repeat(20);
        let mut s = [flank.as_slice(), &ga, &flank, b"C", &b"T".repeat(32), b"C", &flank].concat();
        let n = mask_tandem(&mut s);
        assert_eq!(n, 40 + 32, "{}", String::from_utf8_lossy(&s));
        assert_eq!(&s[..200], flank.as_slice());
        assert!(s[200..240].iter().all(|&b| b == b'N'));
        assert_eq!(&s[240..440], flank.as_slice());
        // runs shorter than a recruitment k-mer stay
        let mut t = [flank.as_slice(), b"C", &b"GA".repeat(10), b"C", &flank].concat();
        assert_eq!(mask_tandem(&mut t), 0);
        // period 6
        let mut u = [flank.as_slice(), &b"ACGTTG".repeat(5), &flank].concat();
        assert!((30..=36).contains(&mask_tandem(&mut u)));
    }

    #[test]
    fn a_microsatellite_in_the_bait_does_not_bring_in_nuclear_microsatellite_reads() {
        // (GA)n in the mitochondrion's seed gap, so only round 2's bait of reads holds it, and 30
        // nuclear reads that are mostly (GA)n
        let dir = tmpdir("ga");
        let mut s = sample(&dir, 0.05);
        let ga = b"GA".repeat(30);
        s.mito.splice(12_000..12_000, ga.iter().copied());
        let mut rng = Rng(45);
        let mut reads = long_reads(&s.mito, 8000, 130, 0.05, "m", &mut rng);
        reads.extend(long_reads(&s.plastid, 8000, 53, 0.05, "p", &mut rng));
        let nuclear = seq(400_000, 33);
        reads.extend((0..30).map(|i| {
            let r = [&nuclear[i * 9000..i * 9000 + 1000], b"GA".repeat(3500).as_slice()].concat();
            (format!("g{i}"), noisy(&r, 0.02, &mut rng))
        }));
        fasta(&s.long, &reads);
        let out = dir.join("out");
        let mt = format!("mitochondrion={}", s.mt_seed.display());
        let pt = format!("plastid={}", s.pt_seed.display());
        let o = options(&["--organelle", "mitochondrion", "--out", out.to_str().unwrap(), "--seeds", &mt, "--seeds", &pt]);
        let m = run(&[s.long.clone()], None, &o).unwrap();
        let recruited = fs::read_to_string(out.join("recruit/round2/mitochondrion.fastq")).unwrap();
        assert!(!recruited.lines().any(|l| l.starts_with("@g")), "nuclear (GA)n reads recruited");
        assert!(m["long_read"]["corrections"][0]["tandem_masked_bases"].as_u64().unwrap() > 40);
        assert_genome(&out, &s.mito, "mitochondrion with a microsatellite");
        fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn deep_reads_are_corrected_down_to_the_depth_cap() {
        let dir = tmpdir("cap");
        let g = seq(5_000, 7);
        let circ = [g.as_slice(), &g[..2000]].concat();
        let deep: Vec<(String, Vec<u8>)> =
            (0..5000).step_by(5).map(|s| (format!("r{s}"), circ[s..s + 2000].to_vec())).collect();
        let path = dir.join("deep.fa");
        fasta(&path, &deep);
        let (reads, cap) = cap_depth(&[path.clone()], &dir.join("capped.fa"), 0).unwrap();
        assert_eq!(reads, vec![dir.join("capped.fa")]);
        let depth = cap["depth"].as_f64().unwrap();
        assert!((350.0..=400.0).contains(&depth), "{cap}");
        let n = cap["kept_reads"].as_u64().unwrap() as f64;
        let expected = 1000.0 * CORRECT_MAX_DEPTH / depth;
        assert!((n - expected).abs() < 0.15 * expected, "{n} reads kept, {expected} expected");
        assert_eq!(read_seqs(&reads[0]).unwrap().len() as f64, n);
        // shallow reads go in as they are
        fasta(&path, &deep.iter().step_by(20).cloned().collect::<Vec<_>>());
        let (reads, cap) = cap_depth(&[path.clone()], &dir.join("capped2.fa"), 0).unwrap();
        assert_eq!(reads, vec![path]);
        assert_eq!(cap["kept_fraction"], json!(1.0));
        assert!(!dir.join("capped2.fa").exists());
        fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn an_organelle_at_the_recruitment_cap_takes_no_further_round() {
        // the plastid's 3 Mb of reads pass a 2.5 Mb cap in round 1; the mitochondrion's seed gap
        // still needs round 2
        let dir = tmpdir("capped");
        let s = sample(&dir, 0.05);
        let out = dir.join("out");
        let mt = format!("mitochondrion={}", s.mt_seed.display());
        let pt = format!("plastid={}", s.pt_seed.display());
        let o = options(&[
            "--organelle", "both", "--out", out.to_str().unwrap(), "--seeds", &mt, "--seeds", &pt,
            "--recruit-max-bases", "2500000",
        ]);
        let report = run(&[s.long.clone()], None, &o).unwrap();
        assert_eq!(report["failed"], json!([]));
        let p = read_json(&out.join("plastid/summary.json")).unwrap();
        let rounds = p["long_read"]["rounds"].as_array().unwrap();
        assert_eq!(rounds[0]["sampled_by"], "max_bases", "{rounds:?}");
        assert_eq!(rounds[0]["used"], json!(true));
        assert_eq!(rounds[1]["used"], json!(false));
        assert_eq!(rounds[1]["capped_before"], json!(true));
        assert_eq!(rounds[1]["bait"], "seeds + corrected reads", "competitive read bait for a capped organelle");
        let purposes: Vec<&Value> =
            p["long_read"]["corrections"].as_array().unwrap().iter().map(|c| &c["purpose"]).collect();
        assert_eq!(purposes, [&json!("assembly")]);
        assert!(p["long_read"]["assembled_reads"].as_str().unwrap().ends_with("correct/plastid.fasta"));
        assert_genome(&out.join("plastid"), &s.plastid, "capped plastid");
        let m = read_json(&out.join("mitochondrion/summary.json")).unwrap();
        assert_eq!(m["long_read"]["rounds"][1]["used"], json!(true));
        assert_genome(&out.join("mitochondrion"), &s.mito, "mitochondrion beside a capped plastid");
        fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn a_round_is_used_only_when_it_adds_a_twentieth() {
        assert!(round_grew(None, 0));
        assert!(round_grew(Some(1000), 1050));
        assert!(!round_grew(Some(1000), 1049));
        assert!(!round_grew(Some(1000), 900));
        assert!(!round_grew(Some(0), 0), "nothing twice");
        assert!(round_grew(Some(0), 1));
    }

    #[test]
    fn a_bait_correction_failure_does_not_discard_the_other_organelle() {
        let dir = tmpdir("bait-isolation");
        let s = sample(&dir, 0.05);
        let out = dir.join("out");
        // A directory at the output FASTA path forces an I/O failure only for MT.
        fs::create_dir_all(out.join("correct/round1/mitochondrion.bait.fasta")).unwrap();
        let mt = format!("mitochondrion={}", s.mt_seed.display());
        let pt = format!("plastid={}", s.pt_seed.display());
        let o = options(&["--organelle", "both", "--out", out.to_str().unwrap(), "--seeds", &mt, "--seeds", &pt]);
        let err = run(&[s.long.clone()], None, &o).unwrap_err();
        assert!(format!("{err:#}").contains("mitochondrion failed"));
        let summary = read_json(&out.join("summary.json")).unwrap();
        assert_eq!(summary["failed"], json!(["mitochondrion"]));
        assert!(summary["organelles"]["mitochondrion"]["error"].as_str().unwrap().contains("bait correction"));
        assert!(out.join("plastid/summary.json").is_file());
        fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn assembly_correction_failure_is_isolated_with_and_without_short_reads() {
        for hybrid in [false, true] {
            let dir = tmpdir(if hybrid { "hybrid-isolation" } else { "self-isolation" });
            let s = sample(&dir, 0.05);
            let out = dir.join("out");
            let blocked = if hybrid { "correct/hybrid/mitochondrion.fasta" } else { "correct/mitochondrion.fasta" };
            fs::create_dir_all(out.join(blocked)).unwrap();
            let mt = format!("mitochondrion={}", s.mt_seed.display());
            let pt = format!("plastid={}", s.pt_seed.display());
            let o = options(&["--organelle", "both", "--out", out.to_str().unwrap(), "--seeds", &mt, "--seeds", &pt]);
            let (s1, s2) = (dir.join("short_1.fa"), dir.join("short_2.fa"));
            let mut a = Vec::new(); let mut b = Vec::new();
            for (tag, genome) in [("m", &s.mito), ("p", &s.plastid)] {
                let circ = [genome.as_slice(), &genome[..400]].concat();
                for (i, st) in (0..genome.len()).step_by(20).enumerate() {
                    let f = &circ[st..st+400];
                    a.push((format!("{tag}{i}"), f[..150].to_vec()));
                    b.push((format!("{tag}{i}"), crate::evidence::revcomp(&f[250..])));
                }
            }
            fasta(&s1, &a); fasta(&s2, &b);
            let short_files = [s1, s2];
            let short = hybrid.then_some(Short { reads: &short_files, max_bases: 0 });
            let err = run(&[s.long.clone()], short, &o).unwrap_err();
            assert!(format!("{err:#}").contains("mitochondrion failed"));
            let summary = read_json(&out.join("summary.json")).unwrap();
            assert_eq!(summary["failed"], json!(["mitochondrion"]));
            let reason = summary["organelles"]["mitochondrion"]["error"].as_str().unwrap();
            assert!(reason.contains(if hybrid { "hybrid correction" } else { "self correction for assembly" }), "{reason}");
            assert!(out.join("plastid/summary.json").is_file());
            fs::remove_dir_all(&dir).ok();
        }
    }

    #[test]
    fn missing_and_nonfile_recruit_outputs_are_errors_but_empty_files_are_valid() {
        let dir = tmpdir("required_recruit_outputs");
        let p = dir.join("numt_candidates.fastq");
        let seeds = BTreeMap::new();
        assert!(read_seqs(&p).is_err());
        assert!(split_partial(&dir, 21, &seeds, &[], 0, 1).is_err());
        fs::create_dir(&p).unwrap();
        assert!(read_seqs(&p).is_err());
        assert!(split_partial(&dir, 21, &seeds, &[], 0, 1).is_err());
        fs::remove_dir(&p).unwrap();
        fs::write(&p, b"").unwrap();
        assert!(read_seqs(&p).unwrap().is_empty());
        assert!(split_partial(&dir, 21, &seeds, &[], 0, 1).unwrap().is_empty());
        fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn hybrid_masking_failure_is_isolated() {
        let dir = tmpdir("hybrid_mask_failure");
        let s = sample(&dir, 0.05);
        let out = dir.join("out");
        fs::create_dir_all(out.join("correct/mitochondrion.masked.fasta")).unwrap();
        let mt = format!("mitochondrion={}", s.mt_seed.display());
        let pt = format!("plastid={}", s.pt_seed.display());
        let o = options(&["--organelle", "both", "--out", out.to_str().unwrap(), "--seeds", &mt, "--seeds", &pt]);
        let s1 = dir.join("short.fa");
        let short_reads = vec![("s0".to_string(), s.plastid[..150].to_vec())];
        fasta(&s1, &short_reads);
        let short_files = [s1];
        let err = run(&[s.long.clone()], Some(Short { reads: &short_files, max_bases: 0 }), &o).unwrap_err();
        assert!(format!("{err:#}").contains("mitochondrion failed"));
        let report = read_json(&out.join("summary.json")).unwrap();
        assert_eq!(report["failed"], json!(["mitochondrion"]));
        assert!(report["organelles"]["mitochondrion"]["error"].as_str().unwrap().contains("masking corrected reads"));
        assert!(out.join("plastid/summary.json").is_file());
        fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn noisy_long_reads_give_both_organelles_and_the_second_round_closes_the_seed_gap() {
        let dir = tmpdir("ont");
        let s = sample(&dir, 0.05);
        let out = dir.join("out");
        let mt = format!("mitochondrion={}", s.mt_seed.display());
        let pt = format!("plastid={}", s.pt_seed.display());
        let o = options(&["--organelle", "both", "--out", out.to_str().unwrap(), "--seeds", &mt, "--seeds", &pt]);
        let report = run(&[s.long.clone()], None, &o).unwrap();
        assert_eq!(report["failed"], json!([]));
        assert_genome(&out.join("mitochondrion"), &s.mito, "mitochondrion");
        assert_genome(&out.join("plastid"), &s.plastid, "plastid");
        let m = read_json(&out.join("mitochondrion/summary.json")).unwrap();
        assert_eq!(m["read_type"], "ont");
        assert_eq!(m["read_set"], "whole_genome");
        assert_eq!(m["seed_database"]["source"], "custom");
        let rounds = m["long_read"]["rounds"].as_array().unwrap();
        assert_eq!(rounds.len(), 2);
        assert_eq!(rounds[1]["used"], json!(true), "the reads across the gap came in round 2");
        assert!(rounds[1]["recruited_bases"].as_u64() > rounds[0]["recruited_bases"].as_u64());
        assert!(m["long_read"]["assembled_reads"].as_str().unwrap().ends_with("correct/mitochondrion.fasta"));
        // nothing nuclear was recruited: a read every 130 bp of a 30 kb genome
        assert!(m["recruited_reads"].as_u64().unwrap() <= 231, "{}", m["recruited_reads"]);
        assert_eq!(m["short_read_evidence"], Value::Null);
        assert!(out.join("mitochondrion/assessment.json").is_file());
        assert!(out.join("summary.json").is_file());
        fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn a_seed_of_genes_only_covers_no_read_and_the_partly_covered_reads_take_over() {
        // the mitochondrial seed is five 2.5 kb "genes" every 6 kb: no 8 kb read is 80% seed
        // (a multi-species library and a plant mitochondrion; ZM4's ONT: none in 12 Gb)
        let dir = tmpdir("genes");
        let s = sample(&dir, 0.05);
        let genes: Vec<(String, Vec<u8>)> = (0..5)
            .map(|i| (format!("gene{i}"), s.mito[i * 6000..i * 6000 + 2500].to_vec()))
            .collect();
        fasta(&s.mt_seed, &genes);
        let out = dir.join("out");
        let mt = format!("mitochondrion={}", s.mt_seed.display());
        let pt = format!("plastid={}", s.pt_seed.display());
        let o = options(&["--organelle", "mitochondrion", "--out", out.to_str().unwrap(), "--seeds", &mt, "--seeds", &pt]);
        let m = run(&[s.long.clone()], None, &o).unwrap();
        let rounds = m["long_read"]["rounds"].as_array().unwrap();
        assert_eq!(rounds[0]["recruited_reads"], json!(0), "{rounds:?}");
        assert!(rounds[1]["recruited_reads"].as_u64().unwrap() > 200, "{rounds:?}");
        assert_eq!(m["long_read"]["corrections"][0]["purpose"], "bait");
        assert!(m["long_read"]["corrections"][0]["partial_reads"].as_u64().unwrap() > 200);
        assert_genome(&out, &s.mito, "mitochondrion from genes");
        // with one round there is nothing to assemble, and the run says why
        let o = options(&[
            "--organelle", "mitochondrion", "--out", dir.join("out1").to_str().unwrap(), "--seeds", &mt,
            "--seeds", &pt, "--recruit-rounds", "1",
        ]);
        let e = run(&[s.long.clone()], None, &o).unwrap_err().to_string();
        assert!(e.contains("more --recruit-rounds"), "{e}");
        fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn one_round_leaves_the_seed_gap_open() {
        let dir = tmpdir("ont1");
        let s = sample(&dir, 0.05);
        let out = dir.join("out");
        let mt = format!("mitochondrion={}", s.mt_seed.display());
        let pt = format!("plastid={}", s.pt_seed.display());
        let o = options(&[
            "--organelle", "mitochondrion", "--out", out.to_str().unwrap(), "--seeds", &mt,
            "--seeds", &pt, "--recruit-rounds", "1",
        ]);
        let m = run(&[s.long.clone()], None, &o).unwrap();
        assert_eq!(m["long_read"]["rounds"].as_array().unwrap().len(), 1);
        assert_ne!(m["circular"], json!([true]), "no read crosses 11,600-12,400: {m}");
        fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn hybrid_reads_correct_with_the_short_reads_and_check_the_graph_with_both() {
        let dir = tmpdir("hybrid");
        let s = sample(&dir, 0.08);
        // 2 x 150 bp pairs of 400 bp fragments, every 20 bp of each genome; nuclear ones every 400
        let mut m1 = Vec::new();
        let mut m2 = Vec::new();
        for (tag, g, step) in [("m", &s.mito, 20), ("p", &s.plastid, 20)] {
            let circ = [g.as_slice(), &g[..400]].concat();
            for (i, st) in (0..g.len()).step_by(step).enumerate() {
                let f = &circ[st..st + 400];
                m1.push((format!("{tag}{i}"), f[..150].to_vec()));
                m2.push((format!("{tag}{i}"), crate::evidence::revcomp(&f[250..])));
            }
        }
        let (s1, s2) = (dir.join("short_1.fa"), dir.join("short_2.fa"));
        fasta(&s1, &m1);
        fasta(&s2, &m2);
        let out = dir.join("out");
        let mt = format!("mitochondrion={}", s.mt_seed.display());
        let pt = format!("plastid={}", s.pt_seed.display());
        let o = options(&["--organelle", "mitochondrion", "--out", out.to_str().unwrap(), "--seeds", &mt, "--seeds", &pt]);
        let short = Short { reads: &[s1, s2], max_bases: 0 };
        run(&[s.long.clone()], Some(short), &o).unwrap();
        assert_genome(&out, &s.mito, "hybrid mitochondrion");
        let m = read_json(&out.join("summary.json")).unwrap();
        assert_eq!(m["read_type"], "ont_hybrid");
        let corrections = m["long_read"]["corrections"].as_array().unwrap();
        assert_eq!(corrections.last().unwrap()["hybrid"], json!(true));
        // by path components: the separator is the platform's
        let assembled = Path::new(m["long_read"]["assembled_reads"].as_str().unwrap());
        assert!(assembled.ends_with(Path::new("hybrid").join("mitochondrion.fasta")));
        assert!(m["long_read"]["short_reads"]["kept_reads"].as_u64().unwrap() >= 2 * 1400);
        // the homopolymer runs are judged from the short reads, not the uncorrected long ones
        assert_eq!(m["homopolymer_polish"]["reads_used"], "short reads");
        assert!(out.join("evidence.short.json").is_file());
        assert!(m["short_read_evidence"]["links"].as_u64().is_some());
        fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn target_reads_of_both_organelles_and_zero_rounds_are_refused() {
        let dir = tmpdir("refuse");
        let reads = dir.join("r.fa");
        fasta(&reads, &[("r".into(), seq(1000, 1))]);
        let out = dir.join("out");
        let o = options(&["--organelle", "both", "--out", out.to_str().unwrap(), "--read-set", "target"]);
        assert!(run(&[reads.clone()], None, &o).unwrap_err().to_string().contains("one organelle"));
        let o = options(&["--organelle", "plastid", "--out", out.to_str().unwrap(), "--recruit-rounds", "0"]);
        assert!(run(&[reads], None, &o).unwrap_err().to_string().contains("at least 1"));
        fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn cross_organelle_partial_filter_drops_downsampled_donor_leaks_and_rescues_crossing_reads() {
        // Mito (12 kb) and plastid (16 kb) share a 5 kb MTPT; plastid is 20x deeper in the WGS
        // library (600x vs 30x) but downsampled to sampling_fraction = 0.2 (120x raw in file).
        let mut mito = seq(12_000, 901);
        let mut plastid = seq(16_000, 902);
        let mtpt = seq(5_000, 903);
        mito[4_000..9_000].copy_from_slice(&mtpt);
        plastid[5_000..10_000].copy_from_slice(&mtpt);
        for p in (4_060..8_940).step_by(120) {
            mito[p] = if mito[p] == b'A' { b'C' } else { b'A' };
        }
        let tile_noisy = |g: &[u8], step: usize, len: usize, seed0: u64| -> Vec<Vec<u8>> {
            let circ = [g, &g[..len]].concat();
            (0..g.len())
                .step_by(step)
                .enumerate()
                .map(|(i, st)| noisy(&circ[st..st + len], 0.05, &mut Rng(seed0 + i as u64)))
                .collect()
        };
        let mt_rec = tile_noisy(&mito, 100, 3_000, 1_000); // ~30x
        let pt_rec = tile_noisy(&plastid, 25, 3_000, 5_000); // ~120x raw, fraction 0.2 => ~600x true

        // 1) Pure plastid read spanning plastid + MTPT (7 kb), misassigned to mito (0): must drop.
        let pure_pt = noisy(&plastid[3_500..10_500], 0.05, &mut Rng(9_001));
        // 2) Genuine crossing mito read: 2.8 kb mt_own + 4.5 kb MTPT (median k-mer sits in MTPT!): must rescue!
        let crossing_mt = noisy(&mito[1_200..8_500], 0.05, &mut Rng(9_002));
        // 3) Pure mt_own read (4 kb): must keep.
        let pure_mt = noisy(&mito[0..4_000], 0.05, &mut Rng(9_003));
        // 4) NUMT + MTPT read: 2.8 kb nuclear + 4.5 kb MTPT: must drop!
        let mut numt_mtpt = seq(2_800, 999);
        numt_mtpt.extend_from_slice(&mito[4_000..8_500]);
        let numt_mtpt = noisy(&numt_mtpt, 0.05, &mut Rng(9_004));

        let partial = [
            (pure_pt.as_slice(), 0usize),
            (crossing_mt.as_slice(), 0usize),
            (pure_mt.as_slice(), 0usize),
            (numt_mtpt.as_slice(), 0usize),
        ];
        let keep = filter_cross_organelle_partials(
            &[mt_rec, pt_rec],
            &[1.0, 0.2],
            &partial,
            2_000,
        );
        assert_eq!(keep, vec![false, true, true, false], "keep={keep:?}");
    }

    #[test]
    fn ont_r10_fastq_auto_runs_phased_overlap_ec_while_hybrid_keeps_legacy() {
        let dir = tmpdir("ont_r10_e2e");
        // 24 kb circular mitochondrion carrying two 1.8 kb diverged repeat copies (6 SNPs apart)
        // and 3 forward-strand CpG methylation hotspots, sequenced as R10.4.1 FASTQ (~1.2% noise).
        let mut mito = seq(24_000, 3_101);
        let copy1 = seq(1_800, 3_102);
        let mut copy2 = copy1.clone();
        for &p in &[200usize, 500, 800, 1_100, 1_400, 1_650] {
            copy2[p] = if copy2[p] == b'A' { b'G' } else { b'A' };
        }
        mito[2_000..3_800].copy_from_slice(&copy1);
        mito[13_000..14_800].copy_from_slice(&copy2);
        for &mp in &[6_000usize, 10_000, 19_000] {
            mito[mp] = b'C';
            mito[mp + 1] = b'G';
        }

        let plastid = seq(16_000, 3_103);
        let mut rng = Rng(3_104);
        let sim_fastq = |genome: &[u8], len: usize, step: usize, tag: &str, meth: &[usize], rng: &mut Rng| -> String {
            let circ = [genome, &genome[..len]].concat();
            let mut buf = String::new();
            for (i, st) in (0..genome.len()).step_by(step).enumerate() {
                let is_fwd = i % 3 != 0;
                let mut truth = circ[st..st + len].to_vec();
                if is_fwd {
                    for &mp in meth {
                        for off in [mp, mp + genome.len()] {
                            if off >= st && off < st + len && rng.unit() < 0.70 {
                                truth[off - st] = b'T';
                            }
                        }
                    }
                }
                let (mut s, mut q) = (Vec::with_capacity(len + 32), Vec::with_capacity(len + 32));
                for &b in &truth {
                    let u = rng.unit();
                    if u < 0.006 {
                        let alt = if b == b'A' { b'C' } else { b'A' };
                        s.push(alt);
                        q.push(b'+'); // Q10
                    } else if u < 0.010 {
                        s.push(b);
                        q.push(b'9'); // Q24
                        s.push(b);
                        q.push(b','); // Q11
                    } else if u < 0.012 {
                        continue;
                    } else {
                        s.push(b);
                        q.push(b'9'); // Q24
                    }
                }
                if !is_fwd {
                    s = crate::evidence::revcomp(&s);
                    q.reverse();
                }
                buf.push_str(&format!(
                    "@{tag}{i}\n{}\n+\n{}\n",
                    std::str::from_utf8(&s).unwrap(),
                    std::str::from_utf8(&q).unwrap()
                ));
            }
            buf
        };

        let fq_path = dir.join("r10.fastq");
        let mut fq = sim_fastq(&mito, 6_000, 100, "m", &[6_000, 10_000, 19_000], &mut rng);
        fq.push_str(&sim_fastq(&plastid, 6_000, 120, "p", &[], &mut rng));
        fs::write(&fq_path, fq).unwrap();

        let (mt_seed, pt_seed) = (dir.join("mt_seed.fa"), dir.join("pt_seed.fa"));
        fasta(
            &mt_seed,
            &[
                ("mt_a".into(), mito[..8_000].to_vec()),
                ("mt_b".into(), mito[11_000..].to_vec()),
            ],
        );
        fasta(&pt_seed, &[("pt".into(), plastid.clone())]);

        let out_ont = dir.join("out_ont");
        let mt_arg = format!("mitochondrion={}", mt_seed.display());
        let pt_arg = format!("plastid={}", pt_seed.display());
        let o = options(&[
            "--organelle",
            "mitochondrion",
            "--out",
            out_ont.to_str().unwrap(),
            "--seeds",
            &mt_arg,
            "--seeds",
            &pt_arg,
        ]);
        run(&[fq_path.clone()], None, &o).unwrap();
        assert_genome(&out_ont, &mito, "R10 mitochondrion");
        let m = read_json(&out_ont.join("summary.json")).unwrap();
        let corrs = m["long_read"]["corrections"].as_array().unwrap();
        assert!(corrs.iter().all(|c| c["chemistry"] == "r10"), "{corrs:?}");
        assert_eq!(m["k"], json!(ASSEMBLE_K_R10));
        assert_eq!(m["homopolymer_polish"]["reads_used"], "r10 phased consensus");
        let mol = &molecules(&out_ont)[0];
        assert_eq!(mol.len(), mito.len());
        assert_eq!(contained(&mito, mol), 1.0, "exact 24 kb circular reconstruction");
        assert_eq!(contained(&copy1, mol), 1.0, "copy1 homogenized");
        assert_eq!(contained(&copy2, mol), 1.0, "copy2 homogenized");

        // Verify `run hybrid` on the exact same R10 FASTQ preserves `legacy` EC!
        let (s1, s2) = (dir.join("sr_1.fa"), dir.join("sr_2.fa"));
        let mut a = Vec::new();
        let mut b = Vec::new();
        let circ = [mito.as_slice(), &mito[..400]].concat();
        for (i, st) in (0..mito.len()).step_by(20).enumerate() {
            let f = &circ[st..st + 400];
            a.push((format!("s{i}"), f[..150].to_vec()));
            b.push((format!("s{i}"), crate::evidence::revcomp(&f[250..])));
        }
        fasta(&s1, &a);
        fasta(&s2, &b);
        let out_hyb = dir.join("out_hyb");
        let o_hyb = options(&[
            "--organelle",
            "mitochondrion",
            "--out",
            out_hyb.to_str().unwrap(),
            "--seeds",
            &mt_arg,
            "--seeds",
            &pt_arg,
        ]);
        let short_arr = [s1, s2];
        run(&[fq_path], Some(Short { reads: &short_arr, max_bases: 0 }), &o_hyb).unwrap();
        let mh = read_json(&out_hyb.join("summary.json")).unwrap();
        let hyb_corrs = mh["long_read"]["corrections"].as_array().unwrap();
        assert!(hyb_corrs.iter().all(|c| c["chemistry"] == "legacy"), "{hyb_corrs:?}");

        fs::remove_dir_all(&dir).ok();
    }
}
