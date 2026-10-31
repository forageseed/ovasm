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
use crate::recruit::{classify, peak_rss_mb, unit_hash, Call, ClassifyParams};

/// k of long-read recruitment. `find recruit --preset ont` defaults to 15, which with the
/// built-in seed library (tens of Mb of k-mers, a twentieth of all 15-mers) covers every read
/// by chance: ZM4's ONT library (CRR2374622, R10.4.1) gave 61 Gb of "plastid" reads out of 131 Gb.
/// At 21 a chance hit is a hundred thousand times rarer, and a read of 5% error still keeps a
/// third of its 21-mers.
pub const RECRUIT_K: usize = 21;

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
    /// Without --seeds the built-in land-plant seed library is used, written to <out>/seeds/.
    /// Both organelles are seeded unless --single-seed, so that each claims its own reads.
    #[arg(long = "seeds")]
    seeds: Vec<String>,
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
    #[arg(long, value_parser = clap::value_parser!(u16).range(1..=256))]
    threads: Option<u16>,
}

#[derive(clap::Args)]
pub struct OntArgs {
    /// Long read files (ONT or PacBio CLR; FASTQ/FASTA, optionally .gz); repeat for several.
    #[arg(long = "reads", required = true)]
    reads: Vec<PathBuf>,
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

pub fn run_ont(args: OntArgs) -> Result<()> {
    finish(run(&args.reads, None, &args.options)?, &args.options)
}

pub fn run_hybrid(args: HybridArgs) -> Result<()> {
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
    if !path.is_file() || fs::metadata(path)?.len() == 0 {
        return Ok(out);
    }
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

/// `reads` down to `CORRECT_MAX_DEPTH`: when their median depth (see `read_depths`) is above it,
/// the reads whose id hashes below the share are written to `out` and it is returned in their
/// place. Also the record of it: depth, kept share, reads and bases kept.
fn cap_depth(reads: &[PathBuf], out: &Path, salt: u64) -> Result<(Vec<PathBuf>, Value)> {
    let mut records: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
    for path in reads {
        if !path.is_file() || fs::metadata(path)?.len() == 0 {
            continue;
        }
        let mut reader = needletail::parse_fastx_file(path)
            .with_context(|| format!("cannot read {}", path.display()))?;
        while let Some(record) = reader.next() {
            let record =
                record.with_context(|| format!("malformed record in {}", path.display()))?;
            records.push((record.id().to_vec(), record.seq().into_owned()));
        }
    }
    let depths = read_depths(&records.iter().map(|r| r.1.as_slice()).collect::<Vec<_>>());
    let depth = depth_quantile(&depths, 0.5);
    let Some(d) = depth.filter(|&d| d > CORRECT_MAX_DEPTH) else {
        return Ok((reads.to_vec(), json!({"depth": depth, "kept_fraction": 1.0})));
    };
    let fraction = CORRECT_MAX_DEPTH / d;
    let mut w = BufWriter::new(fs::File::create(out)?);
    let (mut n, mut b) = (0u64, 0u64);
    for (id, seq) in records.iter().filter(|r| unit_hash(&r.0, salt) < fraction) {
        w.write_all(b">")?;
        w.write_all(id)?;
        w.write_all(b"\n")?;
        w.write_all(seq)?;
        w.write_all(b"\n")?;
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
    let mut calls: Vec<(Vec<u8>, Vec<u8>, usize)> = Vec::new();
    let mut bases = vec![0u64; bait.targets.len()];
    if path.is_file() && fs::metadata(&path)?.len() > 0 {
        let mut reader = needletail::parse_fastx_file(&path)
            .with_context(|| format!("cannot read {}", path.display()))?;
        while let Some(record) = reader.next() {
            let record =
                record.with_context(|| format!("malformed record in {}", path.display()))?;
            let seq = record.seq();
            if let (Call::NumtCandidate(t), _) = classify(&seq, &bait, &params) {
                calls.push((record.id().to_vec(), seq.into_owned(), t));
            }
        }
    }
    // only reads at the depth of what their target recruited (see `at_target_depth`)
    for (t, target) in bait.targets.iter().enumerate() {
        let recruited = read_seqs(&dir.join(format!("{}.fastq", target.name)))?;
        let idx: Vec<usize> = (0..calls.len()).filter(|&i| calls[i].2 == t).collect();
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
                calls[i].2 = usize::MAX;
            }
        }
    }
    for c in &calls {
        if c.2 != usize::MAX {
            bases[c.2] += c.1.len() as u64;
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
        for (id, seq, _) in calls.iter().filter(|c| c.2 == t) {
            if fraction < 1.0 && unit_hash(id, salt) >= fraction {
                continue;
            }
            w.write_all(b">")?;
            w.write_all(id)?;
            w.write_all(b"\n")?;
            w.write_all(seq)?;
            w.write_all(b"\n")?;
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
fn correct(reads: &[PathBuf], short: Option<&Path>, fasta: &Path, salt: u64) -> Result<Value> {
    let json = fasta.with_extension("json");
    let (reads, depth_cap) = cap_depth(reads, &fasta.with_extension("capped.fasta"), salt)?;
    let mut argv = Argv::new("assemble correct")
        .each("--reads", &reads)
        .opt("--out", fasta)
        .opt("--out-json", &json);
    if let Some(s) = short {
        argv = argv.opt("--short-reads", s);
    }
    argv.run()?;
    let report = read_json(&json)?;
    let last = report["rounds"].as_array().and_then(|r| r.last()).cloned();
    Ok(json!({
        "reads": path_str(fasta),
        "report": path_str(&json),
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
    }))
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
    let organelles: &[&'static str] = match o.organelle {
        OrganelleChoice::Mitochondrion => &["mitochondrion"],
        OrganelleChoice::Plastid => &["plastid"],
        OrganelleChoice::Both => &["mitochondrion", "plastid"],
    };
    let single = (organelles.len() == 1).then(|| organelles[0]);
    let mut seeds = seed_map(&o.seeds, single)?;
    fs::create_dir_all(&o.out)?;
    let mut builtin = BTreeSet::new();
    if o.read_set == ReadSet::WholeGenome {
        builtin = crate::seeddb::fill_missing(&mut seeds, &o.out.join("seeds"), single, o.single_seed)?;
        if !builtin.is_empty() {
            log(&format!(
                "no --seeds for {}: using the built-in seed library {}",
                builtin.iter().cloned().collect::<Vec<_>>().join(" and "),
                crate::seeddb::VERSION
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
            let mut c = correct(long, None, &fasta, o.salt)?;
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
                    let mut c = correct(&inputs, None, &fasta, o.salt)?;
                    c["purpose"] = "bait".into();
                    c["partial_reads"] = partial.get(organelle).map(|x| x.1).into();
                    let masked = fasta.with_extension("masked.fasta");
                    c["tandem_masked_bases"] = mask_tandem_fasta(&fasta, &masked)?.into();
                    r.bait = Some(masked);
                    r.corrections.push(c);
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
                log(&format!("self-correcting the {organelle} long reads to assemble"));
                let fasta = o.out.join(format!("correct/{organelle}.fasta"));
                let mut c = correct(&r.raw, None, &fasta, o.salt)?;
                c["purpose"] = "assembly".into();
                r.corrected = fasta;
                r.corrections.push(c);
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
        for (name, r) in &reads {
            if let Ok(r) = r {
                let masked = r.corrected.with_extension("masked.fasta");
                mask_tandem_fasta(&r.corrected, &masked)?;
                argv = argv.opt("--seed-reads", seed_spec(name, &masked));
                baited += 1;
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
                let mut c = correct(&r.raw, Some(&dir.join(format!("{organelle}.fastq"))), &fasta, o.salt)?;
                c["purpose"] = "assembly".into();
                r.corrected = fasta;
                r.corrections.push(c);
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
                assemble(o, organelle, &r, &dir, &seeds, &builtin, files.clone(), pipeline_name, t0)
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

/// The HiFi stages on one organelle's corrected reads, then its summary with the long-read
/// record (and the short-read evidence of `run hybrid`).
#[allow(clippy::too_many_arguments)]
fn assemble(
    o: &LongOptions,
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
    log(&format!("assembling the corrected {organelle} reads"));
    let mut summary = pipeline::run(&pipeline::RunConfig {
        reads: vec![r.corrected.clone()],
        organelle: organelle_choice,
        out: dir.to_path_buf(),
        read_type: ReadType::Hifi,
        read_set: ReadSet::Target,
        seeds: own_seeds,
        builtin_seeds: builtin.clone(),
        recruit: pipeline::RecruitBy::Seeds,
        recruit_depth: None,
        recruit_cache: None,
        auto_recruit_cache: false,
        gated_recruit: false,
        recruit_gate_fold: None,
        recruit_nuclear_depth: None,
        recruit_genome_size: None,
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
    })?;
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
    let seed_info = if o.read_set == ReadSet::Target {
        json!({"source": "none"})
    } else {
        json!({
            "source": if builtin.contains(organelle) { "builtin" } else { "custom" },
            "files": seeds.get(organelle).map(|v| v.iter().map(|p| path_str(p)).collect::<Vec<_>>()),
        })
    };
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
        assert_eq!(rounds[1]["bait"], "seeds", "no bait of reads for a capped organelle");
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
        assert!(m["long_read"]["assembled_reads"].as_str().unwrap().ends_with("hybrid/mitochondrion.fasta"));
        assert!(m["long_read"]["short_reads"]["kept_reads"].as_u64().unwrap() >= 2 * 1400);
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
}
