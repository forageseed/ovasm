//! `ovasm run`: reads in, organelle graph and representative molecules out.
//!
//! The same flow as the desktop runner (`src/organelleverse/desktop/ovasm_runner.py`), in the
//! binary: recruit the organelle's reads (seeds, seed-free discovery, or both), stop when fewer
//! than `MIN_RECRUITED_BASES` were found, assemble at the automatic k, score the graph against
//! the reads, unring it, go down the k ladder until the molecules explain half of the graph, and export the
//! unified graph. Reads that came from depth clusters (`--recruit discover|both`) bring what
//! shares the organelle's cluster, so each assembly graph first loses the components that are
//! not the target organelle (`ovasm components`; the graph the later stages and the ladder's
//! half-of-the-graph rule see is the kept one). The stages are the other subcommands, run in this process with the arguments
//! the runner gives them (so their defaults are the subcommands' own), and each leaves its
//! usual report.
//!
//! `--mode standard` adds downsampled replicates of the recruited reads (see `replicate`): the
//! representative graph stays the full-depth one, and every junction of it gets the share of
//! replicates that reproduce it.

use std::collections::{BTreeMap, BTreeSet};
use std::ffi::OsString;
use std::fmt::Write as _;
use std::fs;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{anyhow, bail, Context, Result};
use clap::Parser;
use rustc_hash::{FxHashMap, FxHashSet};
use serde_json::{json, Map, Value};

use crate::assemble::K_LADDER;
use crate::components::{self, ComponentParams, ComponentsReport};
use crate::evidence::{flip, parse_gfa, Oriented};
use crate::panel::Panel;
use crate::recruit::peak_rss_mb;
use crate::replicate::{self, Reproduced};

/// Fewest recruited bases that can hold an organelle genome at all (the smallest are tens of
/// kb). Luzula sylvatica's mitochondrion (v3 holdout) recruited 4 reads, 4 kb, from distant
/// seeds and failed in unringing with "every segment is a repeat"; ordinary runs recruit
/// millions of bases.
pub const MIN_RECRUITED_BASES: u64 = 100_000;

/// What the result is and is not, written into every summary.
const SCOPE: &str =
    "native assembly and read evidence; no independent assembler certification or polishing";

/// Replicates are drawn at this share of the full depth unless `--replicate-depth` says
/// otherwise: two replicates then share half of their reads, and a larger share would make
/// them the same assembly over again.
pub const REPLICATE_MAX_FRACTION: f64 = 0.5;
/// Deepest automatic replicate (more reads only cost time).
pub const REPLICATE_DEPTH_CAP: f64 = 100.0;
/// Shallowest automatic replicate: the depth from which automatic k gave one circle in every
/// thinned Zou accession (12/12 at 30x, see `assemble::run_auto_k`). So standard mode needs a
/// full depth of 60x.
pub const REPLICATE_MIN_DEPTH: f64 = 30.0;

#[derive(Clone, Copy, PartialEq, Eq, Debug, clap::ValueEnum)]
pub enum OrganelleChoice {
    Mitochondrion,
    Plastid,
    /// Both organelles from one recruitment pass, each in its own directory under --out.
    Both,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug, clap::ValueEnum)]
pub enum ReadType {
    /// PacBio HiFi: sparse (closed-syncmer) graph at the automatic k.
    Hifi,
    /// Illumina short reads: dense graph (every k-mer a node), k from the read length.
    Sr,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug, clap::ValueEnum)]
pub enum ReadSet {
    /// Whole-genome sequencing reads: the organelle's reads are recruited first.
    WholeGenome,
    /// Reads of the target organelle already (e.g. from `ovasm recruit`): assembled as given.
    Target,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug, clap::ValueEnum)]
pub enum RecruitBy {
    /// Seed sequences (--seeds), every seeded organelle competing for its own reads.
    Seeds,
    /// No seeds: depth clusters (`discover`) called by their organelle proteins (`panel`).
    Discover,
    /// Seed-recruited reads plus the depth clusters called as the target organelle.
    Both,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug, clap::ValueEnum)]
pub enum Mode {
    /// Recruit, assemble once (with the k-ladder retry), evidence graph.
    Fast,
    /// Fast, then downsampled replicates: each junction's reproduction rate.
    Standard,
}

pub struct RunConfig {
    pub reads: Vec<PathBuf>,
    pub organelle: OrganelleChoice,
    pub out: PathBuf,
    pub read_type: ReadType,
    pub read_set: ReadSet,
    /// Seed FASTA files per organelle (`mitochondrion`, `plastid`).
    pub seeds: BTreeMap<String, Vec<PathBuf>>,
    pub recruit: RecruitBy,
    /// `ovasm recruit --target-depth` (None: keep every recruited read).
    pub recruit_depth: Option<f64>,
    /// Recruit one organelle without seeding the other.
    pub single_seed: bool,
    pub mode: Mode,
    /// Replicates drawn before the stopping rule is asked.
    pub replicates: usize,
    pub max_replicates: usize,
    pub replicate_depth: Option<f64>,
    /// Salt of the recruitment downsampling; replicate i (from 1) uses `salt + i`.
    pub salt: u64,
    /// Sample name of the unified graph (segment names and PanSN paths).
    pub sample: String,
    /// Protein panel FASTA of the component filter (None: the built-in panel).
    pub component_panel: Option<PathBuf>,
}

fn log(message: &str) {
    eprintln!("[ovasm] run: {message}");
}

/// Arguments of a subcommand, as the command line would give them.
struct Argv(Vec<OsString>);

impl Argv {
    fn new(subcommand: &str) -> Self {
        Argv(vec!["ovasm".into(), subcommand.into()])
    }

    fn arg(mut self, a: impl Into<OsString>) -> Self {
        self.0.push(a.into());
        self
    }

    fn opt(self, flag: &str, value: impl Into<OsString>) -> Self {
        self.arg(flag).arg(value)
    }

    fn each(mut self, flag: &str, values: &[PathBuf]) -> Self {
        for v in values {
            self = self.opt(flag, v);
        }
        self
    }

    /// Run the subcommand in this process. Its thread pool is the one `ovasm run` set up.
    fn run(self) -> Result<()> {
        let name = self.0[1].to_string_lossy().into_owned();
        let cli = crate::Cli::try_parse_from(&self.0)
            .map_err(|e| anyhow!("internal: bad arguments for `ovasm {name}`: {e}"))?;
        crate::dispatch(cli.command).with_context(|| format!("`ovasm {name}` failed"))
    }
}

fn read_json(path: &Path) -> Result<Value> {
    serde_json::from_str(
        &fs::read_to_string(path).with_context(|| format!("cannot read {}", path.display()))?,
    )
    .with_context(|| format!("{} is not JSON", path.display()))
}

fn write_json(path: &Path, value: &Value) -> Result<()> {
    fs::write(path, serde_json::to_string_pretty(value)? + "\n")
        .with_context(|| format!("cannot write {}", path.display()))
}

fn path_str(p: &Path) -> String {
    p.display().to_string()
}

fn file_name(p: &Path) -> String {
    p.file_name()
        .map_or_else(|| path_str(p), |n| n.to_string_lossy().into_owned())
}

// ---------------------------------------------------------------------------------------------
// k ladder
// ---------------------------------------------------------------------------------------------

/// The molecules of an attempt must explain at least this share of the assembly graph for the
/// attempt to count as a result. Carex laevigata's mitochondrion (2.77 Mb, 97.7% of it in the
/// graph): k = 1001 to 351 gave no molecule, k = 251 two molecules of 11,560 and 349 bp that
/// share nothing with the published genome. "A molecule" is not enough.
pub const MIN_EXPLAINED: f64 = 0.5;

/// What one assembly attempt gave.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Outcome {
    pub k: usize,
    pub molecules: usize,
    /// Sum of the molecules' lengths.
    pub molecule_bases: u64,
    /// Total length of the assembly graph's unitigs (`assembly.json`); with the component
    /// filter, of the components it kept (the molecules are the target organelle's, and the
    /// nuclear repeats beside it are not what they should explain).
    pub total_length: u64,
}

impl Outcome {
    /// Molecule bases over graph bases (0 without a molecule or a graph); above 1 when
    /// molecules repeat what the graph holds once.
    pub fn explained(&self) -> f64 {
        if self.molecules == 0 || self.total_length == 0 {
            0.0
        } else {
            self.molecule_bases as f64 / self.total_length as f64
        }
    }

    /// At least one molecule, and the molecules explain at least `MIN_EXPLAINED` of the graph.
    pub fn accepted(&self) -> bool {
        self.molecules >= 1 && self.explained() >= MIN_EXPLAINED
    }
}

/// One assembly attempt, and the files of the best attempt so far while others are tried.
pub trait Attempts {
    /// Assemble at `k` (None: the automatic choice), score and unring.
    fn attempt(&mut self, k: Option<usize>) -> Result<Outcome>;
    /// Set the files of the attempt just made aside as the best so far (replacing any earlier
    /// set-aside files). Called before the next attempt overwrites them.
    fn stash(&mut self) -> Result<()>;
    /// Make the set-aside files the result again.
    fn restore(&mut self) -> Result<()>;
    /// The set-aside files are no longer needed.
    fn drop_stash(&mut self) -> Result<()>;
}

#[derive(Debug, PartialEq)]
pub struct Ladder {
    /// Every attempt, in order.
    pub tried: Vec<Outcome>,
    /// Index in `tried` of the attempt whose files are left as the result.
    pub chosen: usize,
    /// The chosen attempt is accepted (see `Outcome::accepted`). When false the output is not
    /// to be trusted: no k explained half of its graph with molecules.
    pub accepted: bool,
}

impl Ladder {
    pub fn result(&self) -> &Outcome {
        &self.tried[self.chosen]
    }
}

/// The k-ladder retry.
///
/// A large k can leave a gap the rescue cannot fill (Zou AT101/AT105 plastids: k=1001 gave 2
/// unitigs and no molecule, k=701 down to 251 all gave the exact plastome), and the automatic k
/// only looks at depth. So when an attempt is not accepted (no molecule, or molecules that
/// explain under `MIN_EXPLAINED` of the graph), the next smaller k of `K_LADDER` is tried, and
/// the first k that is accepted is the result. When none is, the attempt whose molecules explain
/// the largest share of its graph is the result (the earliest of equals), and `accepted` is
/// false. `retry` is false for dense (short-read) assembly, whose k comes from the read length.
pub fn k_ladder(a: &mut impl Attempts, retry: bool) -> Result<Ladder> {
    let mut tried: Vec<Outcome> = Vec::new();
    let (mut best, mut stashed) = (0usize, false);
    let mut k = None;
    loop {
        let o = a.attempt(k)?;
        let i = tried.len();
        if i == 0 || o.explained() > tried[best].explained() {
            best = i;
        }
        tried.push(o);
        match K_LADDER.iter().copied().find(|&x| x < o.k) {
            Some(next) if !o.accepted() && retry => {
                if best == i {
                    a.stash()?;
                    stashed = true;
                }
                k = Some(next);
            }
            _ => break,
        }
    }
    let last = tried.len() - 1;
    let accepted = tried[last].accepted();
    let chosen = if accepted { last } else { best };
    if stashed {
        if chosen == last {
            a.drop_stash()?;
        } else {
            a.restore()?;
        }
    }
    Ok(Ladder {
        tried,
        chosen,
        accepted,
    })
}

/// Files of one attempt, all in the organelle's output directory.
const ATTEMPT_FILES: [&str; 8] = [
    "assembly.gfa",
    "assembly.json",
    "assembly.filtered.gfa",
    "components.json",
    "evidence.json",
    "linearization.json",
    "molecules.fasta",
    "molecules.partial.fasta",
];
const BEST_ATTEMPT_DIR: &str = ".best_attempt";

/// The component filter for graphs of discovered reads: the organelle it keeps and the
/// protein panel that tells the components' genes.
struct ComponentFilter {
    organelle: String,
    panel: Panel,
}

impl ComponentFilter {
    fn new(organelle: &str, panel: Option<&Path>) -> Result<Self> {
        let k = ComponentParams::default().k;
        let panel = match panel {
            Some(f) => Panel::parse(
                &fs::read_to_string(f).with_context(|| format!("cannot read {}", f.display()))?,
                k,
                &[],
            )?,
            None => Panel::embedded(k, &[])?,
        };
        Ok(Self {
            organelle: organelle.to_string(),
            panel,
        })
    }

    /// `graph` without the components that are not the target organelle, written to `kept`;
    /// the full report to `json`.
    fn keep_target(&self, graph: &Path, kept: &Path, json: &Path) -> Result<ComponentsReport> {
        let report = components::run(
            graph,
            &self.organelle,
            &self.panel,
            &ComponentParams::default(),
            kept,
            json,
        )?;
        if report.removed_components > 0 {
            log(&format!(
                "removed {} of {} graph components ({} bp) that are not the {}; see {}",
                report.removed_components,
                report.components.len(),
                report.removed_length,
                self.organelle,
                file_name(json)
            ));
        }
        Ok(report)
    }
}

struct GraphAttempts<'a> {
    reads: &'a [PathBuf],
    dir: &'a Path,
    read_type: ReadType,
    /// Reads found by depth: each graph goes through the component filter first.
    filter: Option<&'a ComponentFilter>,
}

fn molecules_of(linearization: &Value) -> Vec<&Value> {
    linearization["components"]
        .as_array()
        .into_iter()
        .flatten()
        .flat_map(|c| c["best"]["molecules"].as_array().into_iter().flatten())
        .collect()
}

/// `ovasm evidence` as the runner calls it for this read type.
fn evidence_argv(gfa: &Path, reads: &[PathBuf], out: &Path, read_type: ReadType) -> Argv {
    // k, diagonal tolerance, hits per side, anchor: `_EVIDENCE_PRESETS` of `_ovasm.py`; short
    // reads of 100-150 bp can hold few hits per side and anchor 40 bp
    let (side_hits, anchor) = match read_type {
        ReadType::Hifi => ("10", "500"),
        ReadType::Sr => ("5", "40"),
    };
    Argv::new("evidence")
        .opt("--gfa", gfa)
        .opt("--out", out)
        .opt("--k", "21")
        .opt("--diag-tol", "30")
        .opt("--min-side-hits", side_hits)
        .opt("--min-anchor", anchor)
        .opt("--bootstrap", "1000")
        .each("--reads", reads)
}

impl Attempts for GraphAttempts<'_> {
    fn attempt(&mut self, k: Option<usize>) -> Result<Outcome> {
        let d = self.dir;
        let (graph, fasta) = (d.join("assembly.gfa"), d.join("molecules.fasta"));
        let (filtered, components_json) =
            (d.join("assembly.filtered.gfa"), d.join("components.json"));
        // `linearize` writes the FASTA only when it has a molecule: one left by an earlier run
        // into this directory must not pass for this run's (nor a filter report of one)
        for stale in [&fasta, &filtered, &components_json] {
            if stale.exists() {
                fs::remove_file(stale)?;
            }
        }
        log(&match k {
            Some(k) => format!("assembling the read graph at k={k}"),
            None => "assembling the read graph".into(),
        });
        let mut asm = Argv::new("assemble")
            .opt("--out", &graph)
            .opt("--out-json", d.join("assembly.json"))
            .each("--reads", self.reads);
        if let Some(k) = k {
            asm = asm.opt("--k", k.to_string());
        }
        if self.read_type == ReadType::Sr {
            asm = asm.opt("--s", "0");
        }
        asm.run()?;
        let assembly = read_json(&d.join("assembly.json"))?;
        // reads found by depth: the stages go on with the target's components only
        let (graph, total_length) = match self.filter {
            Some(f) => {
                let kept = f.keep_target(&graph, &filtered, &components_json)?;
                (filtered, kept.kept_length as u64)
            }
            None => (graph, assembly["total_length"].as_u64().unwrap_or(0)),
        };
        log("evaluating read support and repeat pairings");
        evidence_argv(&graph, self.reads, &d.join("evidence.json"), self.read_type).run()?;
        log("generating representative sequences and retaining configuration uncertainty");
        Argv::new("linearize")
            .opt("--gfa", &graph)
            .opt("--evidence", d.join("evidence.json"))
            .opt("--out-fasta", &fasta)
            .opt("--out-json", d.join("linearization.json"))
            .opt("--min-reads", "2")
            .run()?;
        let linearization = read_json(&d.join("linearization.json"))?;
        let molecules = molecules_of(&linearization);
        let o = Outcome {
            k: assembly["k"].as_u64().context("assembly.json has no k")? as usize,
            molecules: molecules.len(),
            molecule_bases: molecules.iter().filter_map(|m| m["length"].as_u64()).sum(),
            total_length,
        };
        if !o.accepted() {
            log(&format!(
                "k={} gave {} molecule(s) of {} bp in a graph of {} bp: under {:.0}% explained",
                o.k,
                o.molecules,
                o.molecule_bases,
                o.total_length,
                MIN_EXPLAINED * 100.0
            ));
        }
        Ok(o)
    }

    fn stash(&mut self) -> Result<()> {
        let stash = self.dir.join(BEST_ATTEMPT_DIR);
        if stash.exists() {
            fs::remove_dir_all(&stash)?;
        }
        fs::create_dir_all(&stash)?;
        for f in ATTEMPT_FILES {
            let from = self.dir.join(f);
            if from.exists() {
                fs::rename(&from, stash.join(f))?;
            }
        }
        Ok(())
    }

    fn restore(&mut self) -> Result<()> {
        log("no k was accepted; keeping the attempt whose molecules explain most of its graph");
        let stash = self.dir.join(BEST_ATTEMPT_DIR);
        for f in ATTEMPT_FILES {
            let to = self.dir.join(f);
            if to.exists() {
                fs::remove_file(&to)?;
            }
            if stash.join(f).exists() {
                fs::rename(stash.join(f), &to)?;
            }
        }
        fs::remove_dir_all(&stash)?;
        Ok(())
    }

    fn drop_stash(&mut self) -> Result<()> {
        fs::remove_dir_all(self.dir.join(BEST_ATTEMPT_DIR))?;
        Ok(())
    }
}

// ---------------------------------------------------------------------------------------------
// recruitment
// ---------------------------------------------------------------------------------------------

/// Reads of one organelle, ready to assemble.
struct Recruited {
    reads: Vec<PathBuf>,
    /// None for target reads (nothing was recruited).
    kept: Option<(u64, u64)>,
    /// `seed_database` of the summary: where the reads came from.
    seed_info: Value,
    /// The reads include depth clusters (`--recruit discover|both`): their graph is run
    /// through the component filter.
    discovered: bool,
}

/// Stop before assembly when recruitment left too little to assemble anything.
fn require_enough(organelle: &str, reads: u64, bases: u64, source: &str) -> Result<()> {
    if bases >= MIN_RECRUITED_BASES {
        return Ok(());
    }
    bail!(
        "ovasm recruited {reads} reads ({bases} bp) for the {organelle} from {source}: too \
         little to assemble anything (fewer than {MIN_RECRUITED_BASES} bp). The seeds are \
         probably too distant from this species; supply a closer reference seed (--seeds), use \
         seed-free discovery (--recruit discover), or give the target reads themselves \
         (--read-set target)."
    )
}

/// Concatenate FASTQ/FASTA files into `dest` as FASTQ, each read once (by its name: the
/// header up to the first blank); returns reads and bases written. Read by read: only the
/// names are held.
fn pool_fastq(sources: &[PathBuf], dest: &Path) -> Result<(u64, u64)> {
    if let Some(parent) = dest.parent() {
        fs::create_dir_all(parent)?;
    }
    let mut seen: FxHashSet<Vec<u8>> = FxHashSet::default();
    let mut w = BufWriter::with_capacity(
        1 << 20,
        fs::File::create(dest).with_context(|| format!("cannot write {}", dest.display()))?,
    );
    let (mut reads, mut bases) = (0u64, 0u64);
    for source in sources {
        let size = fs::metadata(source)
            .with_context(|| format!("cannot read {}", source.display()))?
            .len();
        if size == 0 {
            continue; // nothing recruited: an empty file, which the parser rejects
        }
        let mut parser = needletail::parse_fastx_file(source)
            .with_context(|| format!("cannot read {}", source.display()))?;
        while let Some(record) = parser.next() {
            let record =
                record.with_context(|| format!("malformed record in {}", source.display()))?;
            let id = record.id();
            let name = id
                .split(|b| b.is_ascii_whitespace())
                .next()
                .unwrap_or(id)
                .to_vec();
            if !seen.insert(name) {
                continue;
            }
            let seq = record.seq();
            replicate::write_fastq_record(&mut w, id, &seq, record.qual())?;
            reads += 1;
            bases += seq.len() as u64;
        }
    }
    w.flush()?;
    Ok((reads, bases))
}

/// Depth clusters of the reads, each called by its conserved organelle proteins.
fn called_clusters(cfg: &RunConfig) -> Result<Vec<Value>> {
    log("finding organelle reads by depth, without seeds");
    let dir = cfg.out.join("discover");
    Argv::new("discover")
        .opt("--out", &dir)
        .opt("--salt", cfg.salt.to_string())
        .each("--reads", &cfg.reads)
        .run()?;
    let discovery = read_json(&dir.join("discover.json"))?;
    let clusters = discovery["clusters"]
        .as_array()
        .cloned()
        .unwrap_or_default();
    if clusters.is_empty() {
        return Ok(Vec::new());
    }
    let outputs: Vec<PathBuf> = clusters
        .iter()
        .map(|c| c["output"].as_str().map(PathBuf::from))
        .collect::<Option<_>>()
        .context("discover.json: a cluster has no output")?;
    log("telling the depth clusters apart by their conserved organelle proteins");
    Argv::new("panel")
        .opt("--out-json", dir.join("identify.json"))
        .each("--reads", &outputs)
        .run()?;
    let identification = read_json(&dir.join("identify.json"))?;
    let sets = identification["sets"]
        .as_array()
        .context("identify.json has no sets")?;
    if sets.len() != clusters.len() {
        bail!("panel called {} of {} clusters", sets.len(), clusters.len());
    }
    Ok(clusters
        .iter()
        .zip(sets)
        .map(|(c, s)| {
            json!({
                "name": c["name"],
                "reads": c["kept_reads"],
                "bases": s["bases"],
                "estimated_depth": c["estimated_depth"],
                "call": s["call"],
                "output": c["output"],
            })
        })
        .collect())
}

/// The clusters without their file paths, for the summary.
fn cluster_summary(called: &[Value]) -> Value {
    Value::Array(
        called
            .iter()
            .map(|c| {
                let mut m = c.as_object().cloned().unwrap_or_default();
                m.remove("output");
                Value::Object(m)
            })
            .collect(),
    )
}

fn chosen_outputs(called: &[Value], organelle: &str) -> Vec<PathBuf> {
    called
        .iter()
        .filter(|c| c["call"] == organelle)
        .filter_map(|c| c["output"].as_str().map(PathBuf::from))
        .collect()
}

fn other_organelle(organelle: &str) -> &'static str {
    if organelle == "mitochondrion" {
        "plastid"
    } else {
        "mitochondrion"
    }
}

/// Recruit the reads of every requested organelle from whole-genome reads: one `recruit` pass
/// for all seeded organelles and one `discover` + `panel` for all, whichever are asked for.
/// Per organelle, the reads or why there are none.
fn recruit_all(
    cfg: &RunConfig,
    organelles: &[&'static str],
    files: &mut BTreeMap<String, String>,
) -> Result<BTreeMap<&'static str, Result<Recruited>>> {
    let recruit_dir = cfg.out.join("recruit");
    let by_seeds = cfg.recruit != RecruitBy::Discover;
    let by_depth = cfg.recruit != RecruitBy::Seeds;
    let mut report = None;
    if by_seeds {
        let names: Vec<&String> = cfg.seeds.keys().collect();
        log(&format!(
            "recruiting {} reads with seeds for {:?}",
            organelles.join(" and "),
            names
        ));
        let mut argv = Argv::new("recruit")
            .opt("--out", &recruit_dir)
            .opt(
                "--preset",
                match cfg.read_type {
                    ReadType::Hifi => "hifi",
                    ReadType::Sr => "sr",
                },
            )
            .opt("--salt", cfg.salt.to_string())
            .each("--reads", &cfg.reads);
        for (name, paths) in &cfg.seeds {
            for p in paths {
                let mut spec = OsString::from(format!("{name}="));
                spec.push(p);
                argv = argv.opt("--seed", spec);
            }
        }
        // A multi-species marker database's summed length is not a sample genome size, so the
        // depth is only capped when the caller says the seed is one (--recruit-depth).
        if let Some(d) = cfg.recruit_depth {
            argv = argv.opt("--target-depth", d.to_string());
        }
        argv.run()?;
        files.insert(
            "recruitment".into(),
            path_str(&recruit_dir.join("recruit.json")),
        );
        report = Some(read_json(&recruit_dir.join("recruit.json"))?);
    }
    let called = if by_depth {
        let called = called_clusters(cfg)?;
        files.insert(
            "discovery".into(),
            path_str(&cfg.out.join("discover/discover.json")),
        );
        if !called.is_empty() {
            files.insert(
                "identification".into(),
                path_str(&cfg.out.join("discover/identify.json")),
            );
        }
        Some(called)
    } else {
        None
    };

    let mut out = BTreeMap::new();
    for &organelle in organelles {
        let one = || -> Result<Recruited> {
            let mut seed_info = Map::new();
            let mut source = String::from("discovery");
            let mut seeded: Option<(PathBuf, u64, u64)> = None;
            if let Some(report) = &report {
                let seed_files = |o: &str| -> Vec<String> {
                    cfg.seeds
                        .get(o)
                        .map(|v| v.iter().map(|p| file_name(p)).collect())
                        .unwrap_or_default()
                };
                let own = seed_files(organelle);
                source = own.join(", ");
                seed_info.insert("source".into(), "custom".into());
                seed_info.insert("file".into(), own.first().cloned().into());
                seed_info.insert("files".into(), own.into());
                // The other organelle is seeded too, so that its reads are claimed by their own
                // target. Seeded alone, mitochondrial seeds (whole genomes holding
                // plastid-derived MTPT stretches) extended into the plastid genome, whose reads
                // are 10-30x deeper: the Ajuga and Alisma "mitochondria" of the frozen holdout
                // were their plastid genomes.
                let other = other_organelle(organelle);
                seed_info.insert(
                    "competing_seed".into(),
                    if cfg.seeds.contains_key(other) {
                        json!({"organelle": other, "files": seed_files(other)})
                    } else {
                        Value::Null
                    },
                );
                let index = report["targets"]
                    .as_array()
                    .into_iter()
                    .flatten()
                    .position(|t| t["name"] == organelle)
                    .with_context(|| format!("recruit.json has no target {organelle}"))?;
                let o = &report["outputs"][index];
                seeded = Some((
                    recruit_dir.join(format!("{organelle}.fastq")),
                    o["kept_reads"].as_u64().unwrap_or(0),
                    o["kept_bases"].as_u64().unwrap_or(0),
                ));
            } else {
                seed_info.insert("source".into(), "discover".into());
            }
            let (reads, n_reads, n_bases) = match (&seeded, &called) {
                (Some(s), None) => s.clone(),
                // no seeds at all: organelle reads stand out by depth, and conserved proteins
                // say which depth cluster is which organelle (works far from any reference)
                (None, Some(called)) => {
                    seed_info.insert("clusters".into(), cluster_summary(called));
                    let chosen = chosen_outputs(called, organelle);
                    if chosen.is_empty() {
                        let calls: Vec<String> = called
                            .iter()
                            .map(|c| format!("{}: {}", c["name"], c["call"]))
                            .collect();
                        bail!(
                            "ovasm found no depth cluster whose reads encode the conserved \
                             {organelle} proteins (clusters: [{}]); the organelle may be too \
                             shallow to stand out from the nuclear genome, or its reads may \
                             share a depth cluster with the other organelle's; try seeds or \
                             give the target reads",
                            calls.join(", ")
                        );
                    }
                    let pooled = recruit_dir.join(format!("{organelle}.fastq"));
                    let (n, b) = pool_fastq(&chosen, &pooled)?;
                    log(&format!(
                        "{} cluster(s) hold the {organelle}: {n} reads",
                        chosen.len()
                    ));
                    (pooled, n, b)
                }
                // Seeds give a clean core where a relative is in the database; discovery
                // reaches what seeds cannot (v3 holdout, Luzula plastid: seeds took 6,225 reads
                // and missed a 57 kb stretch the library holds at 400x; its cluster has them
                // all). A cluster called "mixed" holds both organelles' reads and is left to
                // the seeds.
                (Some((seed_reads, seed_n, _)), Some(called)) => {
                    let chosen = chosen_outputs(called, organelle);
                    let union = recruit_dir.join(format!("{organelle}.union.fastq"));
                    let mut sources = vec![seed_reads.clone()];
                    sources.extend(chosen.iter().cloned());
                    let (n, b) = pool_fastq(&sources, &union)?;
                    log(&format!(
                        "seeds {seed_n} reads + {} discovered cluster(s) -> {n} reads for the \
                         {organelle}",
                        chosen.len()
                    ));
                    seed_info.insert(
                        "discovery".into(),
                        json!({
                            "clusters": cluster_summary(called),
                            "seed_reads": seed_n,
                            "union_reads": n,
                        }),
                    );
                    (union, n, b)
                }
                (None, None) => unreachable!("recruitment uses seeds, discovery or both"),
            };
            require_enough(organelle, n_reads, n_bases, &source)?;
            Ok(Recruited {
                reads: vec![reads],
                kept: Some((n_reads, n_bases)),
                seed_info: Value::Object(seed_info),
                discovered: called.is_some(),
            })
        };
        out.insert(organelle, one());
    }
    Ok(out)
}

// ---------------------------------------------------------------------------------------------
// replicates
// ---------------------------------------------------------------------------------------------

/// The depth replicates are drawn at, or why there can be none.
///
/// `depth` is the full depth (read bases over the graph's length). Without a request,
/// replicates take half of it, at most `REPLICATE_DEPTH_CAP`, and need `REPLICATE_MIN_DEPTH`;
/// a requested depth is taken as given. Either way a replicate may hold at most
/// `REPLICATE_MAX_FRACTION` of the reads.
pub fn replicate_depth(depth: f64, requested: Option<f64>) -> Result<f64, String> {
    if !(depth.is_finite() && depth > 0.0) {
        return Err("the full-depth graph is empty, so there is no depth to draw from".into());
    }
    let half = depth * REPLICATE_MAX_FRACTION;
    match requested {
        Some(d) if d > half => Err(format!(
            "replicates at {d:.0}x would each hold more than half of the reads ({depth:.0}x in \
             all), too alike to count as independent"
        )),
        Some(d) => Ok(d),
        None if half < REPLICATE_MIN_DEPTH => Err(format!(
            "depth {depth:.0}x is too low: replicates take half of the reads and need \
             {REPLICATE_MIN_DEPTH:.0}x each (a full depth of {:.0}x)",
            REPLICATE_MIN_DEPTH / REPLICATE_MAX_FRACTION
        )),
        None => Ok(half.min(REPLICATE_DEPTH_CAP)),
    }
}

struct ReplicateGraph {
    info: Value,
    set: BTreeSet<u64>,
    kmers: FxHashSet<u64>,
    /// signature -> (from, to, reads) of the replicate's branching links
    junctions: BTreeMap<u64, (String, String, u64, Vec<u8>)>,
}

fn fell_back(reason: String) -> Value {
    log(&format!("standard mode fell back to fast: {reason}"));
    json!({"status": "fell_back_to_fast", "reason": reason})
}

fn link_reads(evidence: &Value, index: usize) -> u64 {
    evidence["links"][index]["reads"].as_u64().unwrap_or(0)
}

/// Downsampled replicates of `reads` against the representative graph in `dir`; the
/// `replicates` entry of the summary.
fn run_replicates(
    cfg: &RunConfig,
    reads: &[PathBuf],
    dir: &Path,
    assembly: &Value,
    start_k: Option<usize>,
    // the graph the representative molecules and `evidence.json` belong to (the kept
    // components of `assembly.gfa` with the component filter), and the filter itself
    reference_graph: &Path,
    filter: Option<&ComponentFilter>,
) -> Result<Value> {
    if cfg.read_type != ReadType::Hifi {
        return Ok(fell_back(
            "replicates are set up for HiFi reads only (the depth rule is not calibrated for \
             short reads)"
                .into(),
        ));
    }
    if !assembly["corrected_reads"].is_null() {
        return Ok(fell_back(
            "depth is so low that the reads were self-corrected for the full-depth assembly; \
             nothing is left to thin"
                .into(),
        ));
    }
    let bases = assembly["bases"].as_f64().unwrap_or(0.0);
    let length = assembly["total_length"].as_f64().unwrap_or(0.0);
    let depth = if length > 0.0 { bases / length } else { 0.0 };
    let target = match replicate_depth(depth, cfg.replicate_depth) {
        Ok(d) => d,
        Err(reason) => return Ok(fell_back(reason)),
    };
    let fraction = target / depth;
    let rdir = dir.join("replicates");
    if rdir.exists() {
        fs::remove_dir_all(&rdir)?;
    }
    let reference = parse_gfa(reference_graph)?;
    let ref_links = replicate::links(&reference);
    let ref_set = replicate::junction_set(&ref_links);
    let ref_kmers = replicate::spelled_kmers(&reference, &ref_links);
    let ref_evidence = read_json(&dir.join("evidence.json"))?;

    let max = cfg.max_replicates.clamp(1, replicate::MAX_REPLICATES);
    let mut done: Vec<ReplicateGraph> = Vec::new();
    // One k for every replicate. Normally the first replicate's automatic k; but when the
    // representative graph needed a smaller k than the automatic one (the k ladder), the
    // replicates start from that k, or half the reads could pick the k that left the gap again.
    let mut k: Option<usize> = start_k;
    let stop = loop {
        let sets: Vec<BTreeSet<u64>> = done.iter().map(|r| r.set.clone()).collect();
        if let Some(reason) = replicate::stop_reason(&sets, cfg.replicates, max) {
            break reason;
        }
        let i = done.len() + 1;
        let salt = cfg.salt.wrapping_add(i as u64);
        let t0 = Instant::now();
        let d = rdir.join(format!("r{i}"));
        fs::create_dir_all(&d)?;
        let sub = d.join("reads.fastq");
        let (n_reads, n_bases) = replicate::downsample(reads, fraction, salt, &sub)?;
        log(&format!(
            "replicate {i}: {n_reads} reads (~{target:.0}x, salt {salt})"
        ));
        let assembled = d.join("assembly.gfa");
        let sub_reads = std::slice::from_ref(&sub);
        let mut asm = Argv::new("assemble")
            .opt("--out", &assembled)
            .opt("--out-json", d.join("assembly.json"))
            .each("--reads", sub_reads);
        // the same k in every replicate, so that their junction sets differ by what the reads
        // say and not by which repeats k spans
        if let Some(k) = k {
            asm = asm.opt("--k", k.to_string());
        }
        asm.run()?;
        let rep_assembly = read_json(&d.join("assembly.json"))?;
        if !rep_assembly["corrected_reads"].is_null() {
            fs::remove_file(&sub)?;
            return Ok(fell_back(format!(
                "replicates at {target:.0}x are too shallow for any k of the ladder without \
                 self-correcting the reads"
            )));
        }
        let rep_k = rep_assembly["k"]
            .as_u64()
            .context("assembly.json has no k")? as usize;
        k = Some(rep_k);
        // the same graph as the representative one: reads found by depth lose the components
        // that are not the target before their junctions are compared
        let graph = match filter {
            Some(f) => {
                let kept = d.join("assembly.filtered.gfa");
                f.keep_target(&assembled, &kept, &d.join("components.json"))?;
                kept
            }
            None => assembled,
        };
        evidence_argv(&graph, sub_reads, &d.join("evidence.json"), cfg.read_type).run()?;
        // the subsample is the reads whose name hashes under the fraction with this salt: it
        // can be drawn again, so it is not kept
        fs::remove_file(&sub)?;
        let g = parse_gfa(&graph)?;
        let links = replicate::links(&g);
        let evidence = read_json(&d.join("evidence.json"))?;
        let set = replicate::junction_set(&links);
        let mut junctions = BTreeMap::new();
        for j in links.iter().filter(|j| j.branching) {
            junctions.entry(j.signature).or_insert((
                g.label(j.from),
                g.label(j.to),
                link_reads(&evidence, j.index),
                j.span.clone(),
            ));
        }
        let previous = done.last().map(|p: &ReplicateGraph| {
            json!({
                "gained": set.difference(&p.set).count(),
                "lost": p.set.difference(&set).count(),
            })
        });
        done.push(ReplicateGraph {
            info: json!({
                "replicate": i,
                "salt": salt,
                "reads": n_reads,
                "bases": n_bases,
                "k": rep_k,
                "unitigs": rep_assembly["unitigs"],
                "total_length": rep_assembly["total_length"],
                "links": links.len(),
                "junctions": set.len(),
                "links_supported": evidence["links_supported"],
                "junctions_of_representative": ref_set.intersection(&set).count(),
                "change_from_previous": previous,
                "graph": path_str(&graph),
                "evidence": path_str(&d.join("evidence.json")),
                "elapsed_seconds": t0.elapsed().as_secs_f64(),
            }),
            kmers: replicate::spelled_kmers(&g, &links),
            set,
            junctions,
        });
    };
    log(&format!("{} replicates: {stop}", done.len()));

    // every link of the representative graph against every replicate
    let n = done.len();
    let mut rows = Vec::new();
    let mut support: FxHashMap<(Oriented, Oriented), f64> = FxHashMap::default();
    let mut tsv = String::from(
        "from\tto\toverlap\tjunction\tsignature\treads\treplicates\tas_junction\treproduced\tsupport\n",
    );
    for j in &ref_links {
        let status: Vec<Reproduced> = done
            .iter()
            .map(|r| replicate::reproduced(j.signature, &j.span, &r.set, &r.kmers))
            .collect();
        let letters: String = status.iter().map(|s| s.letter()).collect();
        let as_junction = status
            .iter()
            .filter(|&&s| s == Reproduced::Junction)
            .count();
        let present = status.iter().filter(|s| s.present()).count();
        let share = present as f64 / n as f64;
        support.insert(canonical(j.from, j.to), share);
        let (from, to) = (reference.label(j.from), reference.label(j.to));
        let reads = link_reads(&ref_evidence, j.index);
        let sig = replicate::signature_hex(j.signature);
        writeln!(
            tsv,
            "{from}\t{to}\t{}\t{}\t{sig}\t{reads}\t{letters}\t{as_junction}\t{present}\t{share:.3}",
            j.overlap, j.branching
        )?;
        rows.push(json!({
            "from": from,
            "to": to,
            "overlap": j.overlap,
            "junction": j.branching,
            "signature": sig,
            "reads": reads,
            "replicates": letters,
            "as_junction": as_junction,
            "reproduced": present,
            "support": share,
        }));
    }
    fs::write(rdir.join("junctions.tsv"), tsv)?;

    // junctions only replicates hold
    let mut only: BTreeMap<u64, Vec<usize>> = BTreeMap::new();
    for (i, r) in done.iter().enumerate() {
        for sig in r.set.difference(&ref_set) {
            only.entry(*sig).or_default().push(i + 1);
        }
    }
    let replicate_only: Vec<Value> = only
        .iter()
        .map(|(sig, seen)| {
            let (from, to, reads, span) = &done[seen[0] - 1].junctions[sig];
            json!({
                "signature": replicate::signature_hex(*sig),
                "replicates": seen,
                "share": seen.len() as f64 / n as f64,
                "example": {"replicate": seen[0], "from": from, "to": to, "reads": reads},
                // spelled by the representative graph without a branch there (a repeat its k
                // spans, a boundary a few bases away), or not in it at all
                "in_representative": match replicate::spelled(span, &ref_kmers) {
                    Some(true) => "contiguous",
                    Some(false) => "absent",
                    None => "not judged",
                },
            })
        })
        .collect();

    let junction_rows: Vec<&Value> = rows.iter().filter(|r| r["junction"] == true).collect();
    let fully = junction_rows
        .iter()
        .filter(|r| r["reproduced"] == json!(n))
        .count();
    let (tagged, untagged) = tag_unified_links(dir, &cfg.sample, &reference.names, &support, n)?;
    Ok(json!({
        "status": "done",
        "estimated_depth": depth,
        "replicate_depth": target,
        "sampling_fraction": fraction,
        "k": k,
        "salt": cfg.salt,
        "requested": cfg.replicates,
        "max": max,
        "run": n,
        "converged": stop.contains("unchanged"),
        "stop_reason": stop,
        "signature": {
            "flank": replicate::FLANK,
            "kmer": replicate::KMER,
            "hash": "FNV-1a 64 + splitmix64 finaliser over the lexicographically smaller of \
                     the junction span and its reverse complement",
        },
        "junctions_of_representative": junction_rows.len(),
        "junctions_reproduced_by_all": fully,
        "per_replicate": done.iter().map(|r| r.info.clone()).collect::<Vec<_>>(),
        "links": rows,
        "replicate_only_junctions": replicate_only,
        "table": path_str(&rdir.join("junctions.tsv")),
        "unified_links_tagged": tagged,
        "unified_links_untagged": untagged,
    }))
}

fn canonical(a: Oriented, b: Oriented) -> (Oriented, Oriented) {
    let rc = (flip(b), flip(a));
    if (a, b) <= rc {
        (a, b)
    } else {
        rc
    }
}

/// Add replicate support to the links of `organelle.gfa`: `rs:f` (share of replicates that
/// reproduce the link) and `rn:i` (replicates run). A unified link is tagged when both its
/// segments are whole assembly segments (`<sample>.<segment>`) joined by a link of the
/// assembly graph; links inside a segment that unification cut stay untagged. Returns the
/// numbers of links tagged and left.
fn tag_unified_links(
    dir: &Path,
    sample: &str,
    names: &[String],
    support: &FxHashMap<(Oriented, Oriented), f64>,
    replicates: usize,
) -> Result<(usize, usize)> {
    let path = dir.join("organelle.gfa");
    let text = fs::read_to_string(&path)?;
    let index: FxHashMap<&str, usize> = names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i))
        .collect();
    let prefix = format!("{sample}.");
    let end = |name: &str, sign: &str| -> Option<Oriented> {
        Some((*index.get(name.strip_prefix(&prefix)?)?, sign == "+"))
    };
    let (mut tagged, mut untagged) = (0, 0);
    let mut out = String::with_capacity(text.len() + 1024);
    for line in text.lines() {
        out.push_str(line);
        let f: Vec<&str> = line.split('\t').collect();
        if f[0] == "L" && f.len() >= 5 {
            let share = match (end(f[1], f[2]), end(f[3], f[4])) {
                (Some(a), Some(b)) => support.get(&canonical(a, b)),
                _ => None,
            };
            match share {
                Some(s) => {
                    write!(out, "\trs:f:{s:.3}\trn:i:{replicates}")?;
                    tagged += 1;
                }
                None => untagged += 1,
            }
        }
        out.push('\n');
    }
    fs::write(&path, out)?;
    Ok((tagged, untagged))
}

// ---------------------------------------------------------------------------------------------
// one organelle, and the run
// ---------------------------------------------------------------------------------------------

/// Assemble one organelle's reads into `dir` and write its `summary.json`.
fn assemble_organelle(
    cfg: &RunConfig,
    organelle: &str,
    recruited: &Recruited,
    dir: &Path,
    mut files: BTreeMap<String, String>,
) -> Result<Value> {
    let t0 = Instant::now();
    fs::create_dir_all(dir)?;
    let filter = if recruited.discovered {
        log(&format!(
            "keeping the graph components that are the {organelle}"
        ));
        Some(ComponentFilter::new(
            organelle,
            cfg.component_panel.as_deref(),
        )?)
    } else {
        None
    };
    let ladder = k_ladder(
        &mut GraphAttempts {
            reads: &recruited.reads,
            dir,
            read_type: cfg.read_type,
            filter: filter.as_ref(),
        },
        // dense (short-read) assembly takes k from the read length, so it is not retried
        cfg.read_type == ReadType::Hifi,
    )?;
    let assembly = read_json(&dir.join("assembly.json"))?;
    let evidence = read_json(&dir.join("evidence.json"))?;
    let linear = read_json(&dir.join("linearization.json"))?;
    // the graph evidence, molecules and the unified graph are about
    let graph = if filter.is_some() {
        dir.join("assembly.filtered.gfa")
    } else {
        dir.join("assembly.gfa")
    };
    let component_filter = match &filter {
        Some(_) => Some(component_filter_summary(&read_json(
            &dir.join("components.json"),
        )?)),
        None => None,
    };
    log("exporting the unified organelle graph");
    Argv::new("unify")
        .opt("--gfa", &graph)
        .opt("--sample", &cfg.sample)
        .opt("--backend", "ovasm")
        .opt("--organelle", organelle)
        .opt("--max-alternatives", "3")
        .opt("--out", dir.join("organelle.gfa"))
        .opt("--out-json", dir.join("organelle.json"))
        .opt("--evidence", dir.join("evidence.json"))
        .opt("--linearization", dir.join("linearization.json"))
        .run()?;
    let replicates = match cfg.mode {
        Mode::Fast => Value::Null,
        Mode::Standard => {
            // the result came from below the automatic k: start the replicates there
            let k = ladder.result().k;
            run_replicates(
                cfg,
                &recruited.reads,
                dir,
                &assembly,
                (ladder.chosen > 0).then_some(k),
                &graph,
                filter.as_ref(),
            )?
        }
    };
    let mode = match (cfg.mode, replicates["status"].as_str()) {
        (Mode::Standard, Some("done")) => "standard",
        _ => "fast",
    };

    let molecules = molecules_of(&linear);
    let n_links = evidence["links"].as_array().map_or(0, Vec::len);
    for (key, name) in [
        ("assembly_graph", "assembly.gfa"),
        ("organelle_graph", "organelle.gfa"),
        ("graph_metadata", "organelle.json"),
        ("evidence", "evidence.json"),
        ("linearization", "linearization.json"),
        ("summary", "summary.json"),
    ] {
        files.insert(key.into(), path_str(&dir.join(name)));
    }
    if component_filter.is_some() {
        files.insert("filtered_assembly_graph".into(), path_str(&graph));
        files.insert(
            "component_filter".into(),
            path_str(&dir.join("components.json")),
        );
    }
    // what no cover explains, as partial paths (not molecules: see `ovasm linearize`)
    let partial_paths = linear["partial_paths"].as_array().map_or(0, Vec::len);
    if partial_paths > 0 {
        if let Some(f) = linear["partial_fasta"].as_str() {
            files.insert("partial_fasta".into(), f.to_string());
        }
    }
    if !molecules.is_empty() {
        files.insert(
            "assembly_fasta".into(),
            path_str(&dir.join("molecules.fasta")),
        );
    }
    let report = json!({
        "backend": "ovasm",
        "organelle": organelle,
        "read_set": match cfg.read_set {
            ReadSet::WholeGenome => "whole_genome",
            ReadSet::Target => "target_reads",
        },
        "read_type": match cfg.read_type {
            ReadType::Hifi => "hifi_only",
            ReadType::Sr => "short_read",
        },
        "seed_database": recruited.seed_info,
        "recruited_reads": recruited.kept.map(|k| k.0),
        "recruited_bases": recruited.kept.map(|k| k.1),
        "k": ladder.result().k,
        "k_tried": ladder
            .tried
            .iter()
            .map(|o| json!({
                "k": o.k,
                "molecules": o.molecules,
                "molecule_bases": o.molecule_bases,
                "total_length": o.total_length,
                "explained": o.explained(),
                "accepted": o.accepted(),
            }))
            .collect::<Vec<_>>(),
        // false: no k explained MIN_EXPLAINED of its graph with molecules; the files are those
        // of the attempt that came closest and the molecules are not to be trusted
        "k_accepted": ladder.accepted,
        "min_explained": MIN_EXPLAINED,
        "kept_attempt": ladder.chosen,
        "unitigs": assembly["unitigs"],
        "supported_links": evidence["links_supported"],
        "links": n_links,
        "molecules": molecules.len(),
        "lengths": molecules.iter().map(|m| m["length"].clone()).collect::<Vec<_>>(),
        "circular": molecules.iter().map(|m| m["circular"].clone()).collect::<Vec<_>>(),
        // paths the reads join where no cover of molecules was found: not molecules
        "partial_paths": partial_paths,
        "component_filter": component_filter,
        "decisive": linear["decisive"],
        "unsolved_components": linear["unsolved_components"],
        "unbridged_anchors": linear["unbridged_anchors"],
        "skipped_anchors": linear["skipped_anchors"],
        "search_truncated": linear["search_truncated"],
        "scope": SCOPE,
        "mode_requested": match cfg.mode {
            Mode::Fast => "fast",
            Mode::Standard => "standard",
        },
        "mode": mode,
        "replicates": replicates,
        "sample": cfg.sample,
        "files": files,
        "elapsed_seconds": t0.elapsed().as_secs_f64(),
        "peak_rss_mb": peak_rss_mb(),
    });
    write_json(&dir.join("summary.json"), &report)?;
    let mut message = format!(
        "{organelle}: {} representative molecules; supported links {}/{n_links}; decisive={}",
        molecules.len(),
        evidence["links_supported"],
        linear["decisive"]
    );
    if molecules.is_empty() {
        message += "; no representative sequence, inspect the graph and evidence";
    } else if !ladder.accepted {
        message += &format!(
            "; NOT TRUSTWORTHY: no k gave molecules explaining {:.0}% of the graph (best: {:.0}% at k={})",
            MIN_EXPLAINED * 100.0,
            ladder.result().explained() * 100.0,
            ladder.result().k
        );
    }
    if let Some(f) = &component_filter {
        if f["removed_components"].as_u64().unwrap_or(0) > 0 {
            message += &format!(
                "; {} graph components ({} bp) removed as not the {organelle}",
                f["removed_components"], f["removed_length"]
            );
        }
    }
    log(&message);
    Ok(report)
}

/// The filter's report as it goes into `summary.json`: the counts and, for each removed
/// component, what it was and why (`components.json` has every component).
fn component_filter_summary(report: &Value) -> Value {
    let removed: Vec<Value> = report["components"]
        .as_array()
        .into_iter()
        .flatten()
        .filter(|c| c["kept"] == false)
        .map(|c| {
            json!({
                "segments": c["segments"].as_array().map_or(0, Vec::len),
                "length": c["length"],
                "depth": c["depth"],
                "target_genes": c["target_genes"].as_array().map_or(0, Vec::len),
                "other_genes": c["other_genes"].as_array().map_or(0, Vec::len),
                "reason": c["reason"],
            })
        })
        .collect();
    json!({
        "reference_depth": report["reference_depth"],
        "kept_components": report["kept_components"],
        "kept_length": report["kept_length"],
        "removed_components": report["removed_components"],
        "removed_length": report["removed_length"],
        "removed": removed,
    })
}

fn validate(cfg: &RunConfig) -> Result<()> {
    for r in &cfg.reads {
        if !r.is_file() {
            bail!(
                "reads not found: {} (a FASTQ/FASTA file, optionally .gz)",
                r.display()
            );
        }
    }
    for (name, paths) in &cfg.seeds {
        if name != "mitochondrion" && name != "plastid" {
            bail!("--seeds: organelle must be mitochondrion or plastid, got {name:?}");
        }
        for p in paths {
            if !p.is_file() {
                bail!("seed FASTA not found: {}", p.display());
            }
        }
    }
    if !(1..=replicate::MAX_REPLICATES).contains(&cfg.replicates)
        || !(cfg.replicates..=replicate::MAX_REPLICATES).contains(&cfg.max_replicates)
    {
        bail!(
            "--replicates must be 1..={} and --max-replicates between it and {}",
            replicate::MAX_REPLICATES,
            replicate::MAX_REPLICATES
        );
    }
    if cfg.replicate_depth.is_some_and(|d| d.is_nan() || d <= 0.0) {
        bail!("--replicate-depth must be positive");
    }
    if cfg.read_set == ReadSet::Target {
        if cfg.organelle == OrganelleChoice::Both {
            bail!(
                "--read-set target takes the reads of one organelle; run mitochondrion and \
                 plastid each on their own reads"
            );
        }
        if cfg.recruit != RecruitBy::Seeds {
            bail!("--recruit {:?} needs --read-set whole-genome", cfg.recruit);
        }
        return Ok(());
    }
    if cfg.recruit != RecruitBy::Seeds && cfg.read_type != ReadType::Hifi {
        bail!("seed-free discovery (--recruit discover|both) supports HiFi reads only");
    }
    if cfg.recruit == RecruitBy::Discover {
        return Ok(());
    }
    let wanted: &[&str] = match cfg.organelle {
        OrganelleChoice::Mitochondrion => &["mitochondrion"],
        OrganelleChoice::Plastid => &["plastid"],
        OrganelleChoice::Both => &["mitochondrion", "plastid"],
    };
    for o in wanted {
        if !cfg.seeds.contains_key(*o) {
            bail!(
                "whole-genome reads need seeds for the {o}: give --seeds {o}=<fasta> (`ovasm \
                 run` takes the built-in seed library when --seeds is left out), or use \
                 --recruit discover, or give the target reads with --read-set target"
            );
        }
    }
    if wanted.len() == 1 && !cfg.single_seed {
        let other = other_organelle(wanted[0]);
        if !cfg.seeds.contains_key(other) {
            bail!(
                "give --seeds {other}=<fasta> as well: the other organelle is seeded too, so \
                 that its reads are claimed by their own target (mitochondrial seeds alone \
                 extend into the much deeper plastid genome through shared MTPT stretches). \
                 --single-seed recruits without it"
            );
        }
    }
    Ok(())
}

/// Run the pipeline; returns what was written to `<out>/summary.json`. With `--organelle both`
/// that is `{"organelles": {name: summary or {"error": ...}}}` and each organelle's own
/// summary is in its directory.
pub fn run(cfg: &RunConfig) -> Result<Value> {
    let t0 = Instant::now();
    validate(cfg)?;
    fs::create_dir_all(&cfg.out)?;
    let organelles: &[&'static str] = match cfg.organelle {
        OrganelleChoice::Mitochondrion => &["mitochondrion"],
        OrganelleChoice::Plastid => &["plastid"],
        OrganelleChoice::Both => &["mitochondrion", "plastid"],
    };
    let mut files = BTreeMap::new();
    let mut recruited: BTreeMap<&'static str, Result<Recruited>> = match cfg.read_set {
        ReadSet::Target => {
            log("using the supplied target reads; whole-genome recruitment was not run");
            if !cfg.seeds.is_empty() {
                log("note: --seeds are not used with --read-set target");
            }
            organelles
                .iter()
                .map(|&o| {
                    (
                        o,
                        Ok(Recruited {
                            reads: cfg.reads.clone(),
                            kept: None,
                            seed_info: json!({"source": "none"}),
                            discovered: false,
                        }),
                    )
                })
                .collect()
        }
        ReadSet::WholeGenome => {
            if cfg.recruit == RecruitBy::Discover && !cfg.seeds.is_empty() {
                log("note: --seeds are not used with --recruit discover");
            }
            recruit_all(cfg, organelles, &mut files)?
        }
    };
    if let [organelle] = organelles {
        let r = recruited
            .remove(organelle)
            .expect("one entry per organelle")?;
        return assemble_organelle(cfg, organelle, &r, &cfg.out, files);
    }
    let mut results = Map::new();
    let mut failed = Vec::new();
    for &organelle in organelles {
        let result = recruited
            .remove(organelle)
            .expect("one entry per organelle")
            .and_then(|r| {
                assemble_organelle(cfg, organelle, &r, &cfg.out.join(organelle), files.clone())
            });
        match result {
            Ok(report) => {
                results.insert(organelle.into(), report);
            }
            Err(e) => {
                let message = format!("{e:#}");
                log(&format!("{organelle} failed: {message}"));
                results.insert(organelle.into(), json!({"error": message}));
                failed.push(organelle);
            }
        }
    }
    let report = json!({
        "backend": "ovasm",
        "organelles": results,
        "failed": failed,
        "elapsed_seconds": t0.elapsed().as_secs_f64(),
        "peak_rss_mb": peak_rss_mb(),
    });
    write_json(&cfg.out.join("summary.json"), &report)?;
    if !failed.is_empty() {
        bail!(
            "{} failed (see {}); the other organelle's results are in its directory",
            failed.join(" and "),
            cfg.out.join("summary.json").display()
        );
    }
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evidence::revcomp;

    // ---- the k ladder, with scripted attempts -------------------------------------------

    /// Attempts whose outcome is scripted: `first_k` is the automatic choice, `outcomes` maps
    /// a k to (molecules, molecule bases); every graph is `GRAPH` bases long.
    struct Scripted {
        first_k: usize,
        outcomes: Vec<(usize, usize, u64)>,
        asked: Vec<Option<usize>>,
        events: Vec<&'static str>,
    }

    const GRAPH: u64 = 1000;

    impl Scripted {
        fn new(first_k: usize, outcomes: &[(usize, usize, u64)]) -> Self {
            Scripted {
                first_k,
                outcomes: outcomes.to_vec(),
                asked: Vec::new(),
                events: Vec::new(),
            }
        }
    }

    impl Attempts for Scripted {
        fn attempt(&mut self, k: Option<usize>) -> Result<Outcome> {
            self.asked.push(k);
            self.events.push("attempt");
            let used = k.unwrap_or(self.first_k);
            let (molecules, molecule_bases) = self
                .outcomes
                .iter()
                .find(|o| o.0 == used)
                .map_or((0, 0), |o| (o.1, o.2));
            Ok(Outcome {
                k: used,
                molecules,
                molecule_bases,
                total_length: GRAPH,
            })
        }
        fn stash(&mut self) -> Result<()> {
            self.events.push("stash");
            Ok(())
        }
        fn restore(&mut self) -> Result<()> {
            self.events.push("restore");
            Ok(())
        }
        fn drop_stash(&mut self) -> Result<()> {
            self.events.push("drop");
            Ok(())
        }
    }

    fn ks(l: &Ladder) -> Vec<usize> {
        l.tried.iter().map(|o| o.k).collect()
    }

    #[test]
    fn a_molecule_explains_half_the_graph_or_it_does_not_count() {
        let o = |molecules, molecule_bases, total_length| Outcome {
            k: 701,
            molecules,
            molecule_bases,
            total_length,
        };
        assert!(o(1, 500, 1000).accepted(), "exactly half is enough");
        assert!(!o(1, 499, 1000).accepted());
        assert!(
            o(1, 2000, 1000).accepted(),
            "repeats counted twice: above one"
        );
        assert!(!o(0, 0, 1000).accepted());
        assert!(!o(2, 0, 1000).accepted());
        assert!(!o(1, 500, 0).accepted(), "an empty graph explains nothing");
        assert_eq!(o(0, 700, 1000).explained(), 0.0);
        assert_eq!(o(2, 250, 1000).explained(), 0.25);
    }

    #[test]
    fn a_first_attempt_that_is_accepted_is_not_retried() {
        let mut a = Scripted::new(1001, &[(1001, 1, 1000), (701, 1, 1000)]);
        let l = k_ladder(&mut a, true).unwrap();
        assert_eq!(ks(&l), [1001]);
        assert_eq!((l.chosen, l.accepted), (0, true));
        assert_eq!(
            a.asked,
            [None],
            "the first attempt leaves k to the assembler"
        );
        assert_eq!(a.events, ["attempt"]);
    }

    #[test]
    fn no_molecule_goes_to_the_next_smaller_k_and_stops_at_the_first_accepted() {
        // the Zou AT101 plastid: nothing at 1001, one circle at 701
        let mut a = Scripted::new(1001, &[(701, 1, 900), (501, 1, 900)]);
        let l = k_ladder(&mut a, true).unwrap();
        assert_eq!(ks(&l), [1001, 701]);
        assert_eq!((l.result().k, l.accepted), (701, true));
        assert_eq!(a.asked, [None, Some(701)]);
        // 1001 explained nothing (ratio 0 = the best so far), so it was set aside once
        assert_eq!(a.events, ["attempt", "stash", "attempt", "drop"]);
    }

    #[test]
    fn molecules_that_explain_little_of_the_graph_are_rejected_and_the_ladder_goes_on() {
        // Carex: nothing at 1001..351, two small molecules at 251 (349 + 11,560 bp of a
        // 2.76 Mb graph, here scaled to 1000), the real ones at 201
        let mut a = Scripted::new(1001, &[(251, 2, 40), (201, 1, 950)]);
        let l = k_ladder(&mut a, true).unwrap();
        assert_eq!(ks(&l), [1001, 701, 501, 351, 251, 201]);
        assert_eq!(
            l.tried.iter().map(|o| o.accepted()).collect::<Vec<_>>(),
            [false, false, false, false, false, true]
        );
        assert_eq!((l.result().k, l.accepted), (201, true));
        assert_eq!(
            a.asked[1..],
            [Some(701), Some(501), Some(351), Some(251), Some(201)]
        );
        assert_eq!(a.events.last(), Some(&"drop"));
        assert!(!a.events.contains(&"restore"));
    }

    #[test]
    fn the_ladder_continues_from_the_automatic_k_not_from_the_top() {
        // low depth: automatic k was 351; only 201 works
        let mut a = Scripted::new(351, &[(201, 1, 900)]);
        let l = k_ladder(&mut a, true).unwrap();
        assert_eq!(ks(&l), [351, 251, 201]);
        assert_eq!(l.result().k, 201);
        assert_eq!(a.asked, [None, Some(251), Some(201)]);
    }

    #[test]
    fn when_no_k_is_accepted_the_attempt_that_explains_most_is_kept_and_flagged() {
        // best ratio in the middle (501: 0.4); the first (0.1) and last (0.3) are worse
        let mut a = Scripted::new(
            801,
            &[
                (801, 1, 100),
                (701, 1, 200),
                (501, 2, 400),
                (351, 1, 300),
                (251, 0, 0),
            ],
        );
        let l = k_ladder(&mut a, true).unwrap();
        assert_eq!(ks(&l), [801, 701, 501, 351, 251, 201, 151, 101]);
        assert!(l.tried.iter().all(|o| !o.accepted()));
        assert_eq!((l.chosen, l.result().k, l.accepted), (2, 501, false));
        assert_eq!(a.events.last(), Some(&"restore"));
        // the best so far was set aside each time a better one came: 801, 701, 501
        assert_eq!(a.events.iter().filter(|e| **e == "stash").count(), 3);
        assert!(!a.events.contains(&"drop"));
    }

    #[test]
    fn when_the_last_attempt_is_the_best_unaccepted_one_its_files_stay() {
        let mut a = Scripted::new(
            251,
            &[(251, 1, 100), (201, 1, 200), (151, 1, 300), (101, 1, 400)],
        );
        let l = k_ladder(&mut a, true).unwrap();
        assert_eq!(ks(&l), [251, 201, 151, 101]);
        assert_eq!((l.result().k, l.accepted), (101, false));
        assert_eq!(a.events.last(), Some(&"drop"), "nothing to put back");
    }

    #[test]
    fn equal_explained_shares_keep_the_earliest_attempt() {
        let mut a = Scripted::new(351, &[(351, 1, 200), (251, 1, 200), (201, 1, 200)]);
        let l = k_ladder(&mut a, true).unwrap();
        assert_eq!((l.chosen, l.result().k, l.accepted), (0, 351, false));
        assert_eq!(a.events.last(), Some(&"restore"));
        assert_eq!(a.events.iter().filter(|e| **e == "stash").count(), 1);
    }

    #[test]
    fn nothing_at_all_keeps_the_first_attempt() {
        let mut a = Scripted::new(1001, &[]);
        let l = k_ladder(&mut a, true).unwrap();
        assert_eq!(ks(&l), K_LADDER);
        assert_eq!((l.chosen, l.accepted), (0, false));
        assert_eq!(a.events.last(), Some(&"restore"));
    }

    #[test]
    fn the_smallest_k_has_nothing_to_fall_back_to() {
        let mut a = Scripted::new(101, &[]);
        let l = k_ladder(&mut a, true).unwrap();
        assert_eq!(ks(&l), [101]);
        assert_eq!((l.chosen, l.accepted), (0, false));
        assert_eq!(a.events, ["attempt"]);
    }

    #[test]
    fn dense_assembly_is_not_retried() {
        // short reads: k comes from the read length (here 127), smaller ladder values exist
        let mut a = Scripted::new(127, &[(101, 1, 900)]);
        let l = k_ladder(&mut a, false).unwrap();
        assert_eq!(ks(&l), [127]);
        assert_eq!((l.chosen, l.accepted), (0, false));
        assert_eq!(a.events, ["attempt"]);
    }

    #[test]
    fn a_k_between_ladder_values_continues_below_it() {
        let mut a = Scripted::new(801, &[(501, 1, 900)]);
        let l = k_ladder(&mut a, true).unwrap();
        assert_eq!(ks(&l), [801, 701, 501]);
    }

    // ---- replicate depth rule -------------------------------------------------------------

    #[test]
    fn replicates_take_half_the_depth_and_need_thirty() {
        assert_eq!(replicate_depth(150.0, None), Ok(75.0));
        assert_eq!(replicate_depth(60.0, None), Ok(30.0));
        assert_eq!(replicate_depth(1000.0, None), Ok(REPLICATE_DEPTH_CAP));
        assert!(replicate_depth(59.0, None).unwrap_err().contains("too low"));
        assert!(replicate_depth(0.0, None).unwrap_err().contains("empty"));
        // a requested depth is taken as given, up to half of the reads
        assert_eq!(replicate_depth(59.0, Some(20.0)), Ok(20.0));
        assert_eq!(replicate_depth(400.0, Some(200.0)), Ok(200.0));
        assert!(replicate_depth(100.0, Some(60.0))
            .unwrap_err()
            .contains("more than half"));
    }

    // ---- small recruitment pieces ---------------------------------------------------------

    #[test]
    fn too_few_recruited_bases_stop_the_run_and_say_why() {
        assert!(require_enough("plastid", 900, MIN_RECRUITED_BASES, "seed.fa").is_ok());
        let e = require_enough("mitochondrion", 4, 4000, "seed.fa")
            .unwrap_err()
            .to_string();
        assert!(e.contains("recruited 4 reads (4000 bp) for the mitochondrion from seed.fa"));
        assert!(e.contains("too little to assemble anything"));
        assert!(e.contains("closer reference seed") && e.contains("seed-free discovery"));
    }

    fn tmpdir(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("ovasm-run-{name}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        d
    }

    #[test]
    fn pooling_keeps_each_read_once_by_its_name() {
        let dir = tmpdir("pool");
        let (a, b, empty) = (dir.join("a.fq"), dir.join("b.fq"), dir.join("empty.fq"));
        fs::write(&a, "@r1 first\nACGT\n+\nIIII\n@r2\nAACC\n+\nIIII\n").unwrap();
        fs::write(&b, "@r2 again\nAACC\n+\nIIII\n@r3\nGGTTA\n+\nIIIII\n").unwrap();
        fs::write(&empty, "").unwrap();
        let out = dir.join("recruit/pooled.fq");
        let (reads, bases) = pool_fastq(&[a, empty, b], &out).unwrap();
        assert_eq!((reads, bases), (3, 13));
        assert_eq!(
            fs::read_to_string(&out).unwrap(),
            "@r1 first\nACGT\n+\nIIII\n@r2\nAACC\n+\nIIII\n@r3\nGGTTA\n+\nIIIII\n"
        );
        fs::remove_dir_all(&dir).ok();
    }

    // ---- the whole flow on synthetic reads ------------------------------------------------

    fn seq(n: usize, seed: u64) -> Vec<u8> {
        let mut s = seed;
        (0..n)
            .map(|_| {
                s = s
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACGT"[(s >> 62) as usize]
            })
            .collect()
    }

    /// Reads tiled over a circular genome, alternating strands.
    fn tile(genome: &[u8], len: usize, step: usize, tag: &str) -> Vec<(String, Vec<u8>)> {
        let circ = [genome, &genome[..len]].concat();
        (0..genome.len())
            .step_by(step)
            .enumerate()
            .map(|(i, st)| {
                let r = circ[st..st + len].to_vec();
                (
                    format!("{tag}{i}"),
                    if i % 2 == 0 { r } else { revcomp(&r) },
                )
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

    fn config(reads: PathBuf, organelle: OrganelleChoice, out: PathBuf) -> RunConfig {
        RunConfig {
            reads: vec![reads],
            organelle,
            out,
            read_type: ReadType::Hifi,
            read_set: ReadSet::Target,
            seeds: BTreeMap::new(),
            recruit: RecruitBy::Seeds,
            recruit_depth: None,
            single_seed: false,
            mode: Mode::Fast,
            replicates: 3,
            max_replicates: 5,
            replicate_depth: None,
            salt: 0,
            sample: "desktop".into(),
            component_panel: None,
        }
    }

    fn is_rotation(circle: &[u8], genome: &[u8]) -> bool {
        if circle.len() != genome.len() {
            return false;
        }
        let doubled = [genome, genome].concat();
        let rc = revcomp(circle);
        doubled
            .windows(genome.len())
            .any(|w| w == circle || w == rc.as_slice())
    }

    fn molecule(dir: &Path) -> Vec<u8> {
        let text = fs::read_to_string(dir.join("molecules.fasta")).unwrap();
        text.lines()
            .filter(|l| !l.starts_with('>'))
            .collect::<String>()
            .into_bytes()
    }

    #[test]
    fn whole_genome_reads_give_both_organelles_from_one_recruitment() {
        // a 30 kb mitochondrion at 80x, a 20 kb plastid at 240x, nuclear reads at 1x
        let dir = tmpdir("both");
        let (mito, plastid, nuclear) = (seq(30_000, 101), seq(20_000, 202), seq(400_000, 303));
        let mut reads = tile(&mito, 4000, 50, "m");
        reads.extend(tile(&plastid, 4000, 17, "p"));
        reads.extend(
            (0..nuclear.len() - 4000)
                .step_by(4000)
                .enumerate()
                .map(|(i, st)| (format!("n{i}"), nuclear[st..st + 4000].to_vec())),
        );
        let wgs = dir.join("wgs.fa");
        fasta(&wgs, &reads);
        let (mt_seed, pt_seed) = (dir.join("mt.fa"), dir.join("pt.fa"));
        fasta(&mt_seed, &[("mt".into(), mito.clone())]);
        fasta(&pt_seed, &[("pt".into(), plastid.clone())]);
        let mut cfg = config(wgs, OrganelleChoice::Both, dir.join("out"));
        cfg.read_set = ReadSet::WholeGenome;
        cfg.seeds.insert("mitochondrion".into(), vec![mt_seed]);
        cfg.seeds.insert("plastid".into(), vec![pt_seed]);
        let report = run(&cfg).unwrap();
        assert_eq!(report["failed"], json!([]));
        for (organelle, genome, n_reads) in
            [("mitochondrion", &mito, 600), ("plastid", &plastid, 1177)]
        {
            let d = cfg.out.join(organelle);
            for f in [
                "assembly.gfa",
                "assembly.json",
                "evidence.json",
                "linearization.json",
                "molecules.fasta",
                "organelle.gfa",
                "organelle.json",
                "summary.json",
            ] {
                assert!(d.join(f).is_file(), "{organelle}/{f}");
            }
            let s = read_json(&d.join("summary.json")).unwrap();
            assert_eq!(s["lengths"], report["organelles"][organelle]["lengths"]);
            assert_eq!(s["recruited_reads"], json!(n_reads), "{organelle}");
            assert_eq!(s["molecules"], json!(1));
            assert_eq!(s["lengths"], json!([genome.len()]));
            assert_eq!(s["circular"], json!([true]));
            assert_eq!(s["read_set"], "whole_genome");
            assert_eq!(s["read_type"], "hifi_only");
            assert_eq!(s["seed_database"]["source"], "custom");
            assert_eq!(
                s["seed_database"]["competing_seed"]["organelle"],
                other_organelle(organelle)
            );
            assert_eq!(s["k_tried"].as_array().unwrap().len(), 1);
            assert_eq!(s["k_tried"][0]["k"], s["k"]);
            assert_eq!(s["k_accepted"], json!(true));
            assert_eq!(s["k_tried"][0]["accepted"], json!(true));
            assert!(s["k_tried"][0]["explained"].as_f64().unwrap() >= MIN_EXPLAINED);
            assert_eq!(s["mode"], "fast");
            // seed-recruited reads: no component filter, the graph is the assembly graph
            assert_eq!(s["component_filter"], Value::Null, "{organelle}");
            assert!(!d.join("assembly.filtered.gfa").exists());
            assert!(is_rotation(&molecule(&d), genome), "{organelle}");
            let unified = fs::read_to_string(d.join("organelle.gfa")).unwrap();
            assert!(unified.contains("\nP\tdesktop#1#"), "a PanSN path");
        }
        assert!(cfg.out.join("recruit/recruit.json").is_file());
        fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn too_little_recruited_stops_before_assembly() {
        // the seed matches nothing in the reads
        let dir = tmpdir("empty");
        let wgs = dir.join("wgs.fa");
        fasta(&wgs, &tile(&seq(20_000, 1), 3000, 100, "x"));
        let (mt_seed, pt_seed) = (dir.join("mt.fa"), dir.join("pt.fa"));
        fasta(&mt_seed, &[("mt".into(), seq(5000, 2))]);
        fasta(&pt_seed, &[("pt".into(), seq(5000, 3))]);
        let mut cfg = config(wgs, OrganelleChoice::Mitochondrion, dir.join("out"));
        cfg.read_set = ReadSet::WholeGenome;
        cfg.seeds
            .insert("mitochondrion".into(), vec![mt_seed.clone()]);
        // the other organelle's seed is asked for
        let e = format!("{:#}", run(&cfg).unwrap_err());
        assert!(e.contains("--seeds plastid=<fasta> as well"), "{e}");
        cfg.seeds.insert("plastid".into(), vec![pt_seed]);
        let e = format!("{:#}", run(&cfg).unwrap_err());
        assert!(
            e.contains("recruited 0 reads (0 bp) for the mitochondrion from mt.fa"),
            "{e}"
        );
        assert!(
            !cfg.out.join("assembly.gfa").exists(),
            "nothing was assembled"
        );
        assert!(!cfg.out.join("summary.json").exists());
        fs::remove_dir_all(&dir).ok();
    }

    /// A circular genome with a 3 kb repeat in two copies (direct), unique stretches between.
    fn repeat_genome() -> Vec<u8> {
        let (a, r, b) = (seq(22_000, 11), seq(3000, 12), seq(17_000, 13));
        [&a[..], &r[..], &b[..], &r[..]].concat()
    }

    #[test]
    fn standard_mode_reports_each_junction_s_share_of_replicates() {
        // 45 kb at about 160x: replicates at 80x; 6 kb reads cross the 3 kb repeat
        let dir = tmpdir("standard");
        let genome = repeat_genome();
        let reads = dir.join("target.fa");
        fasta(&reads, &tile(&genome, 6000, 40, "r"));
        let mut cfg = config(reads, OrganelleChoice::Mitochondrion, dir.join("out"));
        cfg.mode = Mode::Standard;
        let s = run(&cfg).unwrap();
        assert_eq!(s["read_set"], "target_reads");
        assert_eq!(s["recruited_reads"], Value::Null);
        assert_eq!(s["molecules"], json!(1));
        assert_eq!(s["lengths"], json!([genome.len()]));
        assert!(is_rotation(&molecule(&cfg.out), &genome));
        assert_eq!(s["mode"], "standard");
        let rep = &s["replicates"];
        assert_eq!(rep["status"], "done", "{rep}");
        assert_eq!(rep["run"], json!(3));
        assert_eq!(rep["converged"], json!(true));
        assert!((rep["sampling_fraction"].as_f64().unwrap() - 0.5).abs() < 1e-9);
        // the repeat is collapsed: two ways in and two ways out, four junctions
        assert_eq!(rep["junctions_of_representative"], json!(4));
        assert_eq!(rep["junctions_reproduced_by_all"], json!(4));
        let links = rep["links"].as_array().unwrap();
        assert_eq!(links.len(), 4);
        for l in links {
            assert_eq!(l["replicates"].as_str().unwrap().len(), 3);
            assert_eq!(l["support"], json!(1.0), "{l}");
            assert!(l["reads"].as_u64().unwrap() > 0);
        }
        assert_eq!(rep["replicate_only_junctions"], json!([]));
        // replicate draws use different salts and are themselves reported
        let per = rep["per_replicate"].as_array().unwrap();
        assert_eq!(
            per.iter()
                .map(|r| r["salt"].as_u64().unwrap())
                .collect::<Vec<_>>(),
            [1, 2, 3]
        );
        assert!(per.iter().all(|r| r["k"] == rep["k"]));
        assert!(cfg.out.join("replicates/junctions.tsv").is_file());
        assert!(!cfg.out.join("replicates/r1/reads.fastq").exists());
        // every link of the unified graph carries the share
        let unified = fs::read_to_string(cfg.out.join("organelle.gfa")).unwrap();
        let tagged = unified
            .lines()
            .filter(|l| l.starts_with("L\t") && l.contains("\trs:f:1.000\trn:i:3"))
            .count();
        assert_eq!(json!(tagged), rep["unified_links_tagged"]);
        assert!(tagged > 0);
        // the same run again gives the same replicate table (explicit salts only)
        let again = run(&cfg).unwrap();
        assert_eq!(again["replicates"]["links"], rep["links"]);
        fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn standard_mode_falls_back_to_fast_when_the_depth_is_too_low_and_says_so() {
        // 20 kb at 40x: half would be 20x, under the 30x a replicate needs
        let dir = tmpdir("shallow");
        let genome = seq(20_000, 77);
        let reads = dir.join("target.fa");
        fasta(&reads, &tile(&genome, 4000, 100, "r"));
        let mut cfg = config(reads, OrganelleChoice::Plastid, dir.join("out"));
        cfg.mode = Mode::Standard;
        let s = run(&cfg).unwrap();
        assert_eq!(s["molecules"], json!(1));
        assert_eq!(s["mode_requested"], "standard");
        assert_eq!(s["mode"], "fast");
        assert_eq!(s["replicates"]["status"], "fell_back_to_fast");
        assert!(s["replicates"]["reason"]
            .as_str()
            .unwrap()
            .contains("too low"));
        assert!(!cfg.out.join("replicates").exists());
        fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn short_reads_take_the_dense_graph_and_standard_mode_says_it_fell_back() {
        // 12 kb at 75x in 150 bp reads: k is the read length minus 23, and is not retried
        let dir = tmpdir("short");
        let genome = seq(12_000, 55);
        let reads = dir.join("target.fa");
        fasta(&reads, &tile(&genome, 150, 2, "r"));
        let mut cfg = config(reads, OrganelleChoice::Mitochondrion, dir.join("out"));
        cfg.read_type = ReadType::Sr;
        cfg.mode = Mode::Standard;
        let s = run(&cfg).unwrap();
        assert_eq!(s["read_type"], "short_read");
        assert_eq!(s["k"], json!(127));
        assert_eq!(s["k_tried"].as_array().unwrap().len(), 1);
        assert_eq!(s["k_tried"][0]["k"], json!(127));
        assert_eq!(s["k_tried"][0]["molecules"], json!(1));
        assert_eq!(s["k_accepted"], json!(true));
        assert_eq!(s["lengths"], json!([genome.len()]));
        assert_eq!(s["circular"], json!([true]));
        assert!(is_rotation(&molecule(&cfg.out), &genome));
        assert_eq!(s["mode"], "fast");
        assert_eq!(s["replicates"]["status"], "fell_back_to_fast");
        assert!(s["replicates"]["reason"].as_str().unwrap().contains("HiFi"));
        fs::remove_dir_all(&dir).ok();
    }

    /// A deterministic protein of `n` residues, and DNA that translates to it plus a stop
    /// (one codon per residue).
    fn gene(seed: u64) -> Vec<u8> {
        const RESIDUES: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";
        const CODONS: [&[u8; 3]; 20] = [
            b"GCT", b"TGT", b"GAT", b"GAA", b"TTT", b"GGT", b"CAT", b"ATT", b"AAA", b"CTT", b"ATG",
            b"AAT", b"CCT", b"CAA", b"CGT", b"TCT", b"ACT", b"GTT", b"TGG", b"TAT",
        ];
        let mut dna = Vec::new();
        for &r in &protein(200, seed) {
            dna.extend(CODONS[RESIDUES.iter().position(|&x| x == r).unwrap()]);
        }
        dna.extend(b"TAA");
        dna
    }

    fn protein(n: usize, seed: u64) -> Vec<u8> {
        let mut s = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
        (0..n)
            .map(|_| {
                s ^= s << 13;
                s ^= s >> 7;
                s ^= s << 17;
                b"ACDEFGHIKLMNPQRSTVWY"[(s % 20) as usize]
            })
            .collect()
    }

    #[test]
    fn graphs_of_discovered_reads_lose_the_components_that_are_not_the_target() {
        // a 25 kb mitochondrion with two panel genes at 80x, beside a 15 kb circle without
        // genes at 25x (what a depth cluster brings along: nuclear repeats)
        let dir = tmpdir("filter");
        let genome = [seq(8000, 1), gene(1), seq(8000, 2), gene(2), seq(8000, 3)].concat();
        let other = seq(15_000, 9);
        let mut reads = tile(&genome, 4000, 50, "m");
        reads.extend(tile(&other, 4000, 160, "n"));
        let fa = dir.join("reads.fa");
        fasta(&fa, &reads);
        let panel = dir.join("panel.fa");
        fs::write(
            &panel,
            format!(
                ">mito|cox1|Sp_one\n{}\n>mito|atp1|Sp_one\n{}\n",
                String::from_utf8(protein(200, 1)).unwrap(),
                String::from_utf8(protein(200, 2)).unwrap()
            ),
        )
        .unwrap();
        let mut cfg = config(fa.clone(), OrganelleChoice::Mitochondrion, dir.join("out"));
        cfg.component_panel = Some(panel);
        let recruited = Recruited {
            reads: vec![fa],
            kept: None,
            seed_info: json!({"source": "discover"}),
            discovered: true,
        };
        let s = assemble_organelle(&cfg, "mitochondrion", &recruited, &cfg.out, BTreeMap::new())
            .unwrap();
        let filter = &s["component_filter"];
        assert_eq!(filter["kept_components"], json!(1), "{filter}");
        assert!(
            filter["removed_components"].as_u64().unwrap() >= 1,
            "{filter}"
        );
        assert_eq!(filter["removed"][0]["reason"], "depth_below_reference");
        assert_eq!(filter["removed"][0]["target_genes"], json!(0));
        // the molecules are the genome, and the ladder judged them against the kept graph
        assert_eq!(s["molecules"], json!(1), "{s}");
        assert_eq!(s["lengths"], json!([genome.len()]));
        assert!(is_rotation(&molecule(&cfg.out), &genome));
        assert_eq!(s["k_tried"][0]["total_length"], filter["kept_length"]);
        assert!(s["k_tried"][0]["explained"].as_f64().unwrap() >= MIN_EXPLAINED);
        assert_eq!(s["k_accepted"], json!(true));
        let assembly = read_json(&cfg.out.join("assembly.json")).unwrap();
        assert!(
            assembly["total_length"].as_u64() > filter["kept_length"].as_u64(),
            "the removed circle is in the assembly graph only"
        );
        // files: the full graph, the kept one, the report; unify and evidence used the kept one
        let n_segments = |f: &str| {
            fs::read_to_string(cfg.out.join(f))
                .unwrap()
                .lines()
                .filter(|l| l.starts_with("S\t"))
                .count()
        };
        assert!(n_segments("assembly.filtered.gfa") < n_segments("assembly.gfa"));
        let components = read_json(&cfg.out.join("components.json")).unwrap();
        assert!(components["reference_component"].is_u64(), "{components}");
        assert!(components["components"].as_array().unwrap().len() >= 2);
        let unified = fs::read_to_string(cfg.out.join("organelle.gfa")).unwrap();
        for c in components["components"].as_array().unwrap() {
            for name in c["segments"].as_array().unwrap() {
                let line = format!("S\tdesktop.{}\t", name.as_str().unwrap());
                assert_eq!(
                    unified.contains(&line),
                    c["kept"] == true,
                    "segment {name} of a component the filter {}",
                    c["reason"]
                );
            }
        }
        assert_eq!(
            s["files"]["filtered_assembly_graph"],
            json!(path_str(&cfg.out.join("assembly.filtered.gfa")))
        );
        assert_eq!(
            s["files"]["component_filter"],
            json!(path_str(&cfg.out.join("components.json")))
        );
        fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn target_reads_of_both_organelles_at_once_are_refused() {
        let dir = tmpdir("refuse");
        let reads = dir.join("target.fa");
        fasta(&reads, &[("r".into(), seq(100, 1))]);
        let cfg = config(reads, OrganelleChoice::Both, dir.join("out"));
        let e = run(&cfg).unwrap_err().to_string();
        assert!(e.contains("reads of one organelle"), "{e}");
        fs::remove_dir_all(&dir).ok();
    }
}
