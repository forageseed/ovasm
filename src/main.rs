//! ovasm: OrganelleVerse organelle assembler. `ovasm <module> <command>`: run (short, hifi, ont,
//! hybrid), find, assemble, check, solve, pan.

mod assemble;
mod assess;
mod bait;
mod chromosomes;
mod close;
mod components;
mod correct;
mod evidence;
mod frontier;
mod gpu;
mod genotype;
mod hpc;
mod identify;
mod junctions;
mod kmer;
mod label;
mod linearize;
mod longread;
mod pan;
mod pan_external;
mod pan_backend;
mod panel;
mod phase;
mod pipeline;
mod primary;
mod read_context;
mod read_model;
mod reconstruct;
mod recruit;
mod replay;
mod replicate;
mod retained;
mod seeddb;
mod seedfree;
mod thread_evidence;
mod uncertain;
mod unify;
mod unzip;

use std::path::PathBuf;

use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(name = "ovasm", version, about = "OrganelleVerse organelle assembler")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// One step: reads in, organelle graph and molecules out, by data type.
    #[command(subcommand)]
    Run(RunCommand),
    /// Find the organelle's reads in sequencing data.
    #[command(subcommand)]
    Find(FindCommand),
    /// Assemble organelle reads into a graph.
    #[command(subcommand)]
    Assemble(AssembleCommand),
    /// Check an assembly graph against the reads.
    #[command(subcommand)]
    Check(CheckCommand),
    /// Resolve sequences and chromosomes from a checked graph.
    #[command(subcommand)]
    Solve(SolveCommand),
    /// Pangenome graphs of many samples, and samples reconstructed against them.
    #[command(subcommand)]
    Pan(PanCommand),
}

#[derive(Subcommand)]
enum RunCommand {
    /// Illumina short reads: find, assemble on a dense graph (every k-mer a node, k from the
    /// read length, no retry), check, solve.
    Short(RunArgs),
    /// PacBio HiFi: find, assemble on a sparse graph at the automatic k, retried at the next
    /// smaller k of 1001, 701, 501, 351, 251, 201, 151, 101 until the molecules add up to at
    /// least half of the graph's length (when no k does, the attempt that comes closest is kept
    /// and summary.json says k_accepted: false), check, solve; --mode standard adds downsampled
    /// replicates.
    Hifi(RunArgs),
    /// Noisy long reads (ONT, PacBio CLR): one recruitment pass with the seeds (k 21, never
    /// extended from raw reads), self-correction, further rounds baited by the corrected reads
    /// (--recruit-rounds), then the HiFi stages on the corrected reads (`run hifi --read-set
    /// target`: sparse graph, k ladder, check, solve, assess).
    Ont(longread::OntArgs),
    /// Long reads with short reads of the same sample (--long ONT/CLR, --short Illumina): the
    /// long reads are recruited and self-corrected as in `run ont`, seed the short reads'
    /// recruitment, are corrected again with those short reads (hybrid correction, then self
    /// rounds), go through the HiFi stages, and the graph is checked against both: evidence.json
    /// (corrected long reads) and evidence.short.json (short reads).
    Hybrid(longread::HybridArgs),
}

#[derive(Subcommand)]
enum FindCommand {
    /// Recruit mitochondrial/plastid reads from whole-genome sequencing data by seeds.
    Recruit(RecruitArgs),
    /// Seed-free discovery: find organelle reads by k-mer depth, split them into depth clusters.
    Discover(DiscoverArgs),
    /// Say which organelle each set of reads (e.g. a `find discover` cluster) comes from, by
    /// conserved mitochondrial and plastid proteins; no reference genome needed.
    Classify(PanelArgs),
    /// Assemble each read cluster and label unitigs mitochondrion or plastid by gene content;
    /// write the labelled unitigs as the sample's own seeds.
    Label(LabelArgs),
}

#[derive(Subcommand)]
enum AssembleCommand {
    /// Self-correct noisy organelle reads (ONT, CLR) against their own solid k-mers.
    Correct(CorrectArgs),
    /// Assemble HiFi organelle reads on a closed-syncmer sparse de Bruijn graph (GFA out).
    Graph(AssembleArgs),
    /// Keep the connected components of an assembly graph that are the target organelle (by
    /// conserved proteins and depth); for graphs of reads found by `find discover`.
    Filter(ComponentsArgs),
    /// Rewrite a backend's assembly graph as OV-GFA v1: blunt links, PanSN paths, common tags.
    Unify(UnifyArgs),
}

#[derive(Subcommand)]
enum CheckCommand {
    /// Read support for every link, branch and repeat pairing of an assembly graph (GFA).
    Evidence(EvidenceArgs),
    /// Count the complete reads that walk entry -> node -> exit through every repeat node of a
    /// graph. Diagnostic only: the graph is never changed.
    Unzip(UnzipArgs),
    /// Evidence an assembly's own run already holds, no other tool: single-copy depth, copy
    /// number against the molecules' walks, shallow sequence, reads behind every junction used,
    /// graph and read-backed complexity, panel genes per molecule, plastid LSC/SSC/IR.
    Assess(AssessArgs),
}

#[derive(clap::Args)]
struct AssessArgs {
    /// The assembly graph the evidence and linearization are about.
    #[arg(long)]
    gfa: PathBuf,
    #[arg(long)]
    evidence: PathBuf,
    #[arg(long)]
    linearization: PathBuf,
    /// chromosomes.json (`solve chromosomes`): repeat calls and the chromosome count.
    #[arg(long)]
    chromosomes: Option<PathBuf>,
    /// The molecules' FASTA: the built-in panel's genes on each molecule.
    #[arg(long)]
    molecules: Option<PathBuf>,
    #[arg(long, value_parser = ["mitochondrion", "plastid"])]
    organelle: String,
    #[arg(long = "out-json")]
    out_json: PathBuf,
}

fn assess(args: AssessArgs) -> Result<()> {
    let read = |p: &PathBuf| -> Result<serde_json::Value> {
        Ok(serde_json::from_str(&std::fs::read_to_string(p).with_context(|| {
            format!("cannot read {}", p.display())
        })?)?)
    };
    let chromosomes = args.chromosomes.as_ref().map(read).transpose()?;
    let panel = panel::Panel::embedded(6, &[])?;
    let report = assess::assess(
        &assess::Inputs {
            gfa: &args.gfa,
            evidence: &read(&args.evidence)?,
            linearization: &read(&args.linearization)?,
            chromosomes: chromosomes.as_ref(),
            molecules_fasta: args.molecules.as_deref(),
            organelle: &args.organelle,
        },
        Some(&panel),
    )?;
    std::fs::write(&args.out_json, serde_json::to_string_pretty(&report)?)
        .with_context(|| format!("cannot write {}", args.out_json.display()))?;
    eprintln!("[ovasm] assess: {}", assess::summary(&report));
    Ok(())
}

#[derive(Subcommand)]
enum SolveCommand {
    /// Representative linear (or circular) sequence through read-supported phased paths.
    Linearize(LinearizeArgs),
    /// Export a depth-constrained primary path and unplaced graph segments per component.
    /// These assembly records are not certificates of physical chromosome count.
    Primary(PrimaryArgs),
    /// How many chromosomes the molecules are, and which are conformations of one (recombination
    /// across a repeat with reads for both pairings): a JSON report; annotates the molecules'
    /// FASTA headers and the organelle GFA's paths, sequences unchanged.
    Chromosomes(ChromosomesArgs),
}

#[derive(Subcommand)]
enum PanCommand {
    /// Import PGGB GFA or minigraph rGFA as a validated reconstruction reference.
    Import(PanImportArgs),
    /// Build a graph using PGGB or minigraph and normalize it for reconstruction.
    Build(PanBuildArgs),
    /// Merge samples' OV-GFA graphs into one pangenome graph, keeping every configuration path.
    Merge(PanArgs),
    /// Reconstruct one WGS sample using all paths of mitochondrial and plastid reference graphs.
    Reconstruct(ReconstructArgs),
    /// Closest reference genomes of an assembly by k-mer containment (organelle and species).
    Identify(IdentifyArgs),
    /// Retain exact MT/PT sharing evidence from an existing both-organelle assembly.
    Sharing(SharingArgs),
}

#[derive(clap::Args)]
struct ChromosomesArgs {
    /// The assembly graph the evidence and linearization are about.
    #[arg(long)]
    gfa: PathBuf,
    #[arg(long)]
    evidence: PathBuf,
    #[arg(long)]
    linearization: PathBuf,
    /// molecules.fasta: its headers get chromosome=, conformation= and share= (in place).
    #[arg(long)]
    molecules: Option<PathBuf>,
    /// organelle.gfa: its representative paths get CH:i: and CF:i: tags (in place).
    #[arg(long = "organelle-gfa")]
    organelle_gfa: Option<PathBuf>,
    /// The organelle's seed FASTA files: a component without their sequence that is far
    /// shallower than the main one is possibly nuclear or contamination.
    #[arg(long = "seeds")]
    seeds: Vec<PathBuf>,
    /// Whole-library k-mer depth of each segment (segment<TAB>depth); without it, the kept
    /// reads' depth (dp). Either way a component at twice or half the main one's depth is
    /// source_pending (counted only in the range's maximum).
    #[arg(long)]
    depths: Option<PathBuf>,
    #[arg(long = "out-json")]
    out_json: PathBuf,
}

fn chromosomes(args: ChromosomesArgs) -> Result<()> {
    let read = |p: &PathBuf| -> Result<serde_json::Value> {
        Ok(serde_json::from_str(&std::fs::read_to_string(p).with_context(|| {
            format!("cannot read {}", p.display())
        })?)?)
    };
    let depths = match &args.depths {
        Some(p) => {
            let mut d = std::collections::BTreeMap::new();
            for line in std::fs::read_to_string(p)?.lines() {
                let mut f = line.split_whitespace();
                if let (Some(s), Some(v)) = (f.next(), f.next().and_then(|v| v.parse::<f64>().ok())) {
                    d.insert(s.to_string(), v);
                }
            }
            Some(d)
        }
        None => None,
    };
    let seeds = if args.seeds.is_empty() {
        None
    } else {
        Some(chromosomes::seed_kmers(&args.seeds)?)
    };
    let report = chromosomes::call(
        &args.gfa,
        &read(&args.evidence)?,
        &read(&args.linearization)?,
        depths.as_ref(),
        seeds.as_ref(),
    )?;
    if let Some(p) = args.molecules.as_ref().filter(|p| p.exists()) {
        chromosomes::annotate_fasta(p, &report)?;
    }
    if let Some(p) = args.organelle_gfa.as_ref().filter(|p| p.exists()) {
        chromosomes::annotate_gfa(p, &report)?;
    }
    std::fs::write(&args.out_json, serde_json::to_string_pretty(&report)? + "\n")?;
    match report.count {
        Some(n) => eprintln!("[ovasm] chromosomes: {n}"),
        None if report.fragmented_components => eprintln!(
            "[ovasm] chromosome count unresolved; represented graph groups: {} to {} ({})",
            report.range.0,
            report.range.1,
            report.reasons.join("; ")
        ),
        None => eprintln!(
            "[ovasm] chromosomes: {} to {} ({})",
            report.range.0,
            report.range.1,
            report.reasons.join("; ")
        ),
    }
    Ok(())
}

#[derive(clap::Args)]
struct SharingArgs {
    #[arg(long)]
    assembly: PathBuf,
    #[arg(long)]
    out: PathBuf,
    /// Exact word length; divergent homology and transfer direction are not inferred.
    #[arg(long, default_value_t = 31)]
    k: usize,
}

#[derive(clap::Args)]
struct PanImportArgs {
    #[arg(long)]
    input: PathBuf,
    /// Optional file containing additional validated GFA P/W walks.
    #[arg(long)]
    extra_paths: Option<PathBuf>,
    #[arg(long)]
    out: PathBuf,
    #[arg(long, default_value = "external")]
    sample: String,
}

#[derive(Clone, Copy, clap::ValueEnum)]
enum PanBackendChoice { Pggb, Minigraph }

#[derive(clap::Args)]
struct PanBuildArgs {
    #[arg(long, value_enum)]
    backend: PanBackendChoice,
    /// One assembly per sample, SAMPLE=FASTA; first is the minigraph backbone.
    #[arg(long = "assembly", required = true)]
    assemblies: Vec<String>,
    #[arg(long)]
    out: PathBuf,
    #[arg(long, default_value = "external")]
    sample: String,
    #[arg(long, default_value_t = 8)]
    threads: usize,
    /// Executable path; otherwise resolve the selected backend on PATH.
    #[arg(long)]
    executable: Option<PathBuf>,
}

fn build_external_pan(args: PanBuildArgs) -> Result<()> {
    let inputs = args.assemblies.iter().map(|s| {
        let (name, path) = s.split_once('=').context("--assembly requires SAMPLE=FASTA")?;
        Ok((name.to_string(), PathBuf::from(path)))
    }).collect::<Result<Vec<_>>>()?;
    pan_backend::run(&pan_backend::Config {
        backend: match args.backend { PanBackendChoice::Pggb => pan_backend::Backend::Pggb, PanBackendChoice::Minigraph => pan_backend::Backend::Minigraph },
        inputs, out: args.out, threads: args.threads, executable: args.executable, sample: args.sample,
    })?;
    Ok(())
}

#[derive(clap::Args)]
struct ReconstructArgs {
    /// Target WGS read files, including both paired FASTQ mates.
    #[arg(
        long = "reads",
        required_unless_present = "target_reads",
        conflicts_with = "target_reads"
    )]
    reads: Vec<PathBuf>,
    /// Already recruited target files as ORGANELLE=PATH; both organelles required.
    /// Explicit target-input mode preserves an earlier recruitment without rescanning WGS.
    #[arg(long = "target-reads", conflicts_with = "target_depth")]
    target_reads: Vec<String>,
    /// Recruited HiFi reads as ORGANELLE=SAMPLE=PATH; repeat for prior samples.
    #[arg(long = "prior-reads")]
    prior_reads: Vec<String>,
    /// Blunt reference pan GFA as mitochondrion=PATH or plastid=PATH; both required.
    #[arg(long = "reference", required = true)]
    references: Vec<String>,
    #[arg(long)]
    sample: String,
    #[arg(long)]
    out: PathBuf,
    /// Optional resource cap: kept bases / stored reference graph bases (each node once).
    /// This basis is neither eight genomes summed nor certified physical genome size.
    #[arg(long)]
    target_depth: Option<f64>,
    #[arg(long, default_value_t = 0)]
    salt: u64,
    #[arg(long)]
    threads: Option<usize>,
    /// Recruit from the reference graphs alone, without first walking the individual's
    /// sequence that no reference path holds on the WGS k-mer counts.
    #[arg(long)]
    no_frontier: bool,
    /// Memory of that walk's k-mer count table, in GiB (rounded down to a power of two).
    #[arg(long, default_value_t = 4.0)]
    frontier_memory_gib: f64,
}

fn reconstruct(args: ReconstructArgs) -> Result<()> {
    if !(args.frontier_memory_gib >= 1.0 / 1024.0) {
        bail!("frontier memory must be at least 1 MiB");
    }
    // two bytes a counter
    let slots_log2 = ((args.frontier_memory_gib * (1u64 << 30) as f64 / 2.0).log2().floor() as u32).min(34);
    let mut references = std::collections::BTreeMap::new();
    for spec in args.references {
        let (o, p) = spec
            .split_once('=')
            .context("reference must be ORGANELLE=GFA")?;
        if !["mitochondrion", "plastid"].contains(&o)
            || references.insert(o.to_string(), PathBuf::from(p)).is_some()
        {
            bail!("give one mitochondrial and one plastid reference GFA");
        }
    }
    let mut target_reads = std::collections::BTreeMap::<String, Vec<PathBuf>>::new();
    for spec in args.target_reads {
        let (o, p) = spec
            .split_once('=')
            .context("target reads must be ORGANELLE=PATH")?;
        if !["mitochondrion", "plastid"].contains(&o) {
            bail!("invalid target organelle {o}");
        }
        target_reads
            .entry(o.to_string())
            .or_default()
            .push(PathBuf::from(p));
    }
    let mut prior_reads = std::collections::BTreeMap::<String, std::collections::BTreeMap<String, Vec<PathBuf>>>::new();
    for spec in args.prior_reads {
        let (organelle, rest) = spec.split_once('=').context("prior reads must be ORGANELLE=SAMPLE=PATH")?;
        let (sample, path) = rest.split_once('=').context("prior reads must be ORGANELLE=SAMPLE=PATH")?;
        if !["mitochondrion", "plastid"].contains(&organelle) || sample.is_empty() || sample == args.sample || path.is_empty() || !sample.chars().all(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '-' | '.')) {
            bail!("invalid prior sample, organelle, path, or target sample used as its own prior");
        }
        prior_reads.entry(organelle.to_string()).or_default().entry(sample.to_string()).or_default().push(PathBuf::from(path));
    }
    if let Some(n) = args.threads {
        if n == 0 {
            bail!("threads must be positive");
        }
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()?;
    }
    reconstruct::run(&reconstruct::Config {
        reads: args.reads,
        target_reads,
        prior_reads,
        references,
        sample: args.sample,
        out: args.out,
        target_depth: args.target_depth,
        salt: args.salt,
        frontier: (!args.no_frontier).then(|| frontier::Params {
            slots_log2,
            ..Default::default()
        }),
    })?;
    Ok(())
}

#[derive(clap::Args)]
struct IdentifyArgs {
    /// Assembled sequences (FASTA), e.g. `molecules.fasta` from `ovasm solve linearize`.
    #[arg(long)]
    query: PathBuf,
    /// Reference set as NAME=PATH (one genome per record); repeat for several sets.
    #[arg(long = "reference", required = true)]
    references: Vec<String>,
    #[arg(long, default_value_t = 21)]
    k: usize,
    /// Output JSON report.
    #[arg(long)]
    out: PathBuf,
}

fn identify(args: IdentifyArgs) -> Result<()> {
    if !(1..=kmer::MAX_K).contains(&args.k) {
        bail!("--k must be between 1 and {}", kmer::MAX_K);
    }
    let references = args
        .references
        .iter()
        .map(|spec| {
            let (name, path) = spec
                .split_once('=')
                .with_context(|| format!("--reference must be NAME=PATH, got {spec:?}"))?;
            Ok((name.to_string(), PathBuf::from(path)))
        })
        .collect::<Result<Vec<_>>>()?;
    let r = identify::run(&args.query, &references, args.k, &args.out)?;
    if let Some(best) = r.references.first() {
        eprintln!(
            "[ovasm] {} bp, {} k-mers; closest {} ({}): containment {:.3}, identity ~{:.4}",
            r.query_length,
            r.query_kmers,
            best.id,
            best.set,
            best.containment,
            best.estimated_identity
        );
    }
    Ok(())
}

#[derive(clap::Args)]
struct RunArgs {
    /// Read files (FASTQ/FASTA, optionally .gz); repeat for several.
    #[arg(long = "reads", required = true)]
    reads: Vec<PathBuf>,
    /// Organelle to assemble; `both` recruits once and assembles each into <out>/<organelle>/.
    #[arg(long, value_enum)]
    organelle: pipeline::OrganelleChoice,
    /// Output directory: assembly.gfa/json, evidence.json, linearization.json, molecules.fasta
    /// (when there is a molecule), organelle.gfa/json (OV-GFA v1), summary.json; recruit/ and
    /// discover/ for whole-genome reads; replicates/ in standard mode.
    #[arg(long)]
    out: PathBuf,
    /// Set by the subcommand: `run hifi` or `run short`.
    #[arg(skip = pipeline::ReadType::Hifi)]
    read_type: pipeline::ReadType,
    /// whole-genome: recruit the organelle's reads first (stops when fewer than 100 kb are
    /// found). target: the reads are the organelle's already and are assembled as given.
    #[arg(long, value_enum, default_value_t = pipeline::ReadSet::WholeGenome)]
    read_set: pipeline::ReadSet,
    /// Seed FASTA as ORGANELLE=PATH (mitochondrion=mt.fa, plastid=pt.fa; repeat for more
    /// files); a bare PATH seeds the --organelle. Optional: an organelle without --seeds uses
    /// the built-in land-plant seed library (seeddb/, compiled in; it includes Arabidopsis
    /// thaliana), written to <out>/seeds/. Both organelles are always seeded unless
    /// --single-seed, so that each claims its own reads.
    #[arg(long = "seeds")]
    seeds: Vec<String>,
    /// How whole-genome reads are recruited: seeds; discover (no seeds: depth clusters called
    /// by their organelle proteins, HiFi only); both (seed reads plus the clusters called as
    /// the target, each read once). With discovery each assembly graph first loses the
    /// components that are not the target organelle (`ovasm assemble filter`; assembly.filtered.gfa,
    /// components.json), and the later stages work on the kept ones.
    #[arg(long, value_enum, default_value_t = pipeline::RecruitBy::Seeds)]
    recruit: pipeline::RecruitBy,
    /// Downsample seed-recruited reads to about this depth (recruited bases / seed length).
    /// Right for a single close reference (the desktop runner uses 150 for a custom seed);
    /// leave unset for a multi-species seed database, whose summed length is no genome size.
    #[arg(long)]
    recruit_depth: Option<f64>,
    /// New scratch directory for lossless recruitment input replay (disk space required).
    /// Without it, HiFi recruitment caches the reads in <out>/recruit/input-cache when that
    /// disk has room (a quarter of a byte per base and the read names), so that passes after the
    /// first, and the depth sketch, read the cache instead of decompressing the input again.
    #[arg(long)]
    recruit_cache: Option<PathBuf>,
    /// Never cache the reads for recruitment (see --recruit-cache).
    #[arg(long, conflicts_with = "recruit_cache")]
    no_recruit_cache: bool,
    /// HiFi: downsample each organelle's recruited reads to about this read depth before
    /// assembly (their median depth in the depth sketch); 0 keeps every read.
    #[arg(long, default_value_t = 250.0)]
    recruit_max_depth: f64,
    /// Keep at most about this many recruited bases per organelle, which bounds the memory of
    /// assembly whatever recruitment found; 0 for no cap. Short reads too: a plastid at a
    /// seventh of a library (WL354: 35 M reads, 5.3 Gb) went to assembly whole. Assembly memory follows the nuclear
    /// sequence among the reads, not their bases: 4.3 Gb of Medicago polymorpha plastid reads
    /// took 2.5 GB, 8.5 Gb of mostly nuclear Carex echinata "plastid" reads 309 GB. A billion
    /// bases still holds the largest mitochondrial genome (11 Mb) at 90x.
    #[arg(long, default_value_t = 1_000_000_000)]
    recruit_max_bases: u64,
    /// HiFi: the depth gate's organelle threshold for seed recruitment, in multiples of the
    /// nuclear peak depth (`find recruit --gate-fold`; default 3). Extension already lowers it
    /// for an organelle shallower than that but deeper than the nucleus.
    #[arg(long, conflicts_with = "no_extend_gate")]
    extend_gate_fold: Option<f64>,
    /// HiFi: the depth gate's nuclear peak, given instead of judged from the reads' median
    /// depths (`find recruit --nuclear-depth`), which stand two copies deep in a repeat-rich or
    /// polyploid genome.
    #[arg(long, conflicts_with_all = ["no_extend_gate", "genome_size"])]
    nuclear_depth: Option<u32>,
    /// HiFi: the nuclear genome's collapsed size for the depth gate (`find recruit
    /// --genome-size`): all subgenomes of an allopolyploid (1C), one set of an autopolyploid;
    /// bases, or with k, m or g (1.26g).
    #[arg(long, value_parser = parse_genome_size, conflicts_with = "no_extend_gate")]
    genome_size: Option<u64>,
    /// HiFi: recruit without the depth gate (`find recruit` without --extend-gate --bootstrap),
    /// as up to 0.2.1. On a species far from the seed library extension can then run into the
    /// nucleus (Carex echinata: 8.5 Gb of "plastid" reads); --recruit-max-bases still bounds what
    /// goes on to assembly.
    #[arg(long)]
    no_extend_gate: bool,
    /// Seed recruitment: most scanning passes (`ovasm find recruit --iterations`). Default: 60 for
    /// HiFi reads, whose gated recruitment grows a large mitochondrion by a few tens of kilobases
    /// a pass (see pipeline::GATED_RECRUIT_PASSES); `find recruit`'s own default otherwise. A target
    /// stops as soon as it stops growing, so the limit only costs time where one is incomplete.
    #[arg(long)]
    recruit_iterations: Option<usize>,
    /// Checked complete original-read nuclear alignments (ovasm.whole-read-context.v1).
    #[arg(long)]
    read_context: Option<PathBuf>,
    /// Recruit one organelle without a seed for the other (not advised: mitochondrial seeds
    /// alone can extend into the plastid genome).
    #[arg(long)]
    single_seed: bool,
    /// fast: assemble once. standard: also assemble downsampled replicates of the recruited
    /// reads and report, per junction of the full-depth graph, the share of replicates that
    /// reproduce it (summary.json, replicates/junctions.tsv, rs:f/rn:i on organelle.gfa
    /// links). Replicates take half of the reads (at most 100x) and need 30x each, so a full
    /// depth of 60x; below that, and for short reads, the run falls back to fast and says so.
    #[arg(long, value_enum, default_value_t = pipeline::Mode::Fast)]
    mode: pipeline::Mode,
    /// Replicates drawn before the stopping rule applies (1-5). Further ones are drawn until
    /// the junction set was the same in three replicates in a row, or --max-replicates.
    #[arg(long, default_value_t = 3)]
    replicates: usize,
    #[arg(long, default_value_t = replicate::MAX_REPLICATES)]
    max_replicates: usize,
    /// Depth of each replicate instead of half the full depth (at most half of it).
    #[arg(long)]
    replicate_depth: Option<f64>,
    /// Salt of every random draw: recruitment downsampling uses it, replicate i uses salt + i.
    #[arg(long, default_value_t = 0)]
    salt: u64,
    /// Sample name in the unified graph (segment names, PanSN paths).
    #[arg(long, default_value = "sample")]
    sample: String,
    /// Protein panel FASTA (headers `class|gene|species`) of the component filter that
    /// graphs of discovered reads (--recruit discover|both) go through; default: the
    /// built-in panel. See `ovasm assemble filter`.
    #[arg(long)]
    component_panel: Option<PathBuf>,
    #[arg(long, value_parser = clap::value_parser!(u16).range(1..=256))]
    threads: Option<u16>,
}

fn run(args: RunArgs) -> Result<()> {
    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n as usize)
            .build_global()?;
    }
    if args.recruit_max_depth.is_nan() || args.recruit_max_depth < 0.0 {
        bail!("--recruit-max-depth must be 0 (no cap) or positive");
    }
    if args.recruit_iterations == Some(0) {
        bail!("--recruit-iterations must be at least 1");
    }
    if args.extend_gate_fold.is_some_and(|f| f.is_nan() || f <= 0.0) {
        bail!("--extend-gate-fold must be positive");
    }
    if args.nuclear_depth == Some(0) {
        bail!("--nuclear-depth must be positive");
    }
    let gate_options = args.extend_gate_fold.is_some()
        || args.nuclear_depth.is_some()
        || args.genome_size.is_some()
        || args.no_extend_gate;
    if gate_options && args.read_type != pipeline::ReadType::Hifi {
        bail!(
            "--extend-gate-fold, --nuclear-depth, --genome-size and --no-extend-gate set the HiFi depth gate: run hifi only (short reads are gated by their k-mer depths)"
        );
    }
    if gate_options
        && (args.read_set == pipeline::ReadSet::Target
            || args.recruit == pipeline::RecruitBy::Discover)
    {
        bail!(
            "--extend-gate-fold, --nuclear-depth, --genome-size and --no-extend-gate set seed recruitment's depth gate: not with --read-set target or --recruit discover"
        );
    }
    let single = match args.organelle {
        pipeline::OrganelleChoice::Mitochondrion => Some("mitochondrion"),
        pipeline::OrganelleChoice::Plastid => Some("plastid"),
        pipeline::OrganelleChoice::Both => None,
    };
    let mut seeds: std::collections::BTreeMap<String, Vec<PathBuf>> = Default::default();
    for spec in &args.seeds {
        let (name, path) = match (spec.split_once('='), single) {
            (Some((name, path)), _) => (name, path),
            (None, Some(organelle)) => (organelle, spec.as_str()),
            (None, None) => {
                bail!("--seeds must be ORGANELLE=PATH with --organelle both, got {spec:?}")
            }
        };
        seeds.entry(name.to_string()).or_default().push(path.into());
    }
    // Whole-genome reads are recruited with seeds unless discovery does it: an organelle
    // without --seeds takes the built-in library (its other organelle too, unless --single-seed).
    let mut builtin_seeds = std::collections::BTreeSet::new();
    if args.read_set == pipeline::ReadSet::WholeGenome
        && args.recruit != pipeline::RecruitBy::Discover
    {
        let dir = args.out.join("seeds");
        builtin_seeds = seeddb::fill_missing(&mut seeds, &dir, single, args.single_seed)?;
        if !builtin_seeds.is_empty() {
            eprintln!(
                "[ovasm] run: no --seeds for {}: using the built-in seed library {} ({})",
                builtin_seeds
                    .iter()
                    .cloned()
                    .collect::<Vec<_>>()
                    .join(" and "),
                seeddb::VERSION,
                dir.display()
            );
        }
    }
    let report = pipeline::run(&pipeline::RunConfig {
        reads: args.reads,
        organelle: args.organelle,
        out: args.out.clone(),
        read_type: args.read_type,
        read_set: args.read_set,
        seeds,
        builtin_seeds,
        recruit: args.recruit,
        recruit_depth: args.recruit_depth,
        recruit_cache: args.recruit_cache,
        auto_recruit_cache: !args.no_recruit_cache,
        // HiFi recruitment in `run` is gated (`find recruit --extend-gate --bootstrap`) unless
        // --no-extend-gate: the 0.2.1 ungated mode ran into the nucleus on species far from the
        // seed library (see docs/release-0.3.0.md).
        gated_recruit: !args.no_extend_gate,
        recruit_gate_fold: args.extend_gate_fold,
        recruit_nuclear_depth: args.nuclear_depth,
        recruit_genome_size: args.genome_size,
        assemble_k: None,
        recruit_max_depth: (args.recruit_max_depth > 0.0).then_some(args.recruit_max_depth),
        recruit_max_bases: (args.recruit_max_bases > 0).then_some(args.recruit_max_bases),
        recruit_iterations: args.recruit_iterations,
        read_context: args.read_context,
        single_seed: args.single_seed,
        mode: args.mode,
        replicates: args.replicates,
        max_replicates: args.max_replicates,
        replicate_depth: args.replicate_depth,
        salt: args.salt,
        sample: args.sample,
        component_panel: args.component_panel,
    })?;
    eprintln!(
        "[ovasm] run: summary -> {}; peak RSS {}",
        args.out.join("summary.json").display(),
        report["peak_rss_mb"]
            .as_f64()
            .map_or("n/a".into(), |m| format!("{m:.0} MB"))
    );
    Ok(())
}

#[derive(clap::Args)]
struct CorrectArgs {
    /// Organelle reads (e.g. from `ovasm find recruit`), FASTQ/FASTA, optionally .gz.
    #[arg(long = "reads", required = true)]
    reads: Vec<PathBuf>,
    /// k of each correction round, in order (repeat the flag); small k first for noisy reads.
    #[arg(long = "k", default_values_t = vec![17, 25])]
    k: Vec<usize>,
    /// Solid-count threshold (default: the valley after the error peak, per round).
    #[arg(long)]
    min_count: Option<u32>,
    /// Allowed relative length difference between a bridging path and the read gap.
    #[arg(long, default_value_t = 0.2)]
    length_tolerance: f64,
    /// Longest read gap (bp) searched for a bridge.
    #[arg(long, default_value_t = 1000)]
    max_gap: usize,
    /// Graph steps searched per gap.
    #[arg(long, default_value_t = 20_000)]
    max_steps: usize,
    /// Keep read ends outside the first and last solid anchor (trimmed by default).
    #[arg(long)]
    keep_ends: bool,
    /// Solid k-mers must reach this fraction of their best sibling's count (0 disables);
    /// removes systematic minority errors such as homopolymer length calls.
    #[arg(long, default_value_t = 0.0)]
    relative: f64,
    /// Auto solid threshold is at least this fraction of the solid k-mers' peak count.
    #[arg(long, default_value_t = 0.3)]
    peak_fraction: f64,
    /// Accurate short reads (Illumina) whose k-mers define the solid set (hybrid correction);
    /// repeat for several files, e.g. both mates.
    #[arg(long = "short-reads")]
    short_reads: Vec<PathBuf>,
    /// With --short-reads: self-correction rounds (the reads' own k-mers) after the hybrid ones,
    /// for stretches the short reads do not cover (repeat the flag) [default: 15 21 25].
    #[arg(long = "self-k")]
    self_k: Vec<usize>,
    /// With --short-reads: hybrid rounds only.
    #[arg(long)]
    no_self_rounds: bool,
    #[arg(long)]
    threads: Option<usize>,
    /// Output FASTA of corrected reads.
    #[arg(long)]
    out: PathBuf,
    #[arg(long)]
    out_json: PathBuf,
}

fn correct(args: CorrectArgs) -> Result<()> {
    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()?;
    }
    let r = correct::run(
        args.reads,
        args.short_reads,
        &args.k,
        &if args.no_self_rounds {
            Vec::new()
        } else if args.self_k.is_empty() {
            vec![15, 21, 25]
        } else {
            args.self_k.clone()
        },
        correct::CorrectParams {
            k: args.k[0],
            min_count: args.min_count,
            length_tolerance: args.length_tolerance,
            max_gap: args.max_gap,
            max_steps: args.max_steps,
            keep_ends: args.keep_ends,
            relative: args.relative,
            peak_fraction: args.peak_fraction,
        },
        &args.out,
        &args.out_json,
    )?;
    eprintln!(
        "[ovasm] corrected reads -> {}; {:.1}s",
        r.output, r.elapsed_seconds
    );
    Ok(())
}

#[derive(clap::Args)]
struct AssembleArgs {
    /// Record cleaning proposals without deleting unitigs (diagnostic, default off).
    #[arg(long)]
    cleaning_audit_only: bool,
    /// Organelle reads (e.g. from `ovasm find recruit`), FASTQ/FASTA, optionally .gz.
    #[arg(long = "reads", required = true)]
    reads: Vec<PathBuf>,
    /// Node k-mer length (odd). Longer k gives cleaner graphs at organelle depth; 1001 and 701
    /// both reach 16/16 Nipponbare junctions (tolerant scoring), 501 only 10-14. Default: the
    /// largest of 1001, 701, ... 101 whose node depth reaches 22 (see `assemble::run_auto_k`);
    /// in dense mode the median read length minus 23 (see `assemble::dense_k`).
    #[arg(long)]
    k: Option<usize>,
    /// Syncmer s-mer length; 0 makes every k-mer a node (dense mode, for short reads).
    #[arg(long, default_value_t = 31)]
    s: usize,
    /// Node depth automatic k must reach (default 22; see `assemble::run_auto_k`).
    #[arg(long)]
    auto_k_depth: Option<f64>,
    /// Minimum reads for a node or an edge (default: a tenth of the typical node depth, >= 3).
    #[arg(long)]
    min_count: Option<u32>,
    /// Automatic minimum reads as this fraction of the typical node depth (default 0.1, at
    /// least 3; in (0, 1]). Ignored when --min-count is given. A lower fraction keeps
    /// low-frequency connections but also more error nodes at high depth; on one HiFi sample
    /// 0.03-0.05 recovered minor repeat-mediated junctions, and on 17 other HiFi samples
    /// 4 of 17 gained low-depth fragments or extra sequence, so it is opt-in.
    #[arg(long)]
    min_count_fraction: Option<f64>,
    /// Nuclear depth over the organelle's own depth in the whole library (`ovasm run` passes
    /// recruitment's): a dead-end branch off a junction no deeper than 1.5 times that share of
    /// the node depth is a NUMT's nuclear flank and is removed. Not given: no such removal.
    #[arg(long)]
    nuclear_share: Option<f64>,
    /// Complete original-read nuclear context; shared local matches do not establish source.
    #[arg(long)]
    read_context: Option<PathBuf>,
    /// Legacy name of `--rust-hpc` kept for existing callers; not a chromosome-count setting.
    #[arg(long, hide = true, conflicts_with = "auto_k_depth")]
    multichrom_graph: bool,
    /// Experimental homopolymer-run coordinates for the same assembly engine: reads are
    /// run-compressed as they are read in, run lengths restored from original reads, then the
    /// shared junction stage runs. --k counts runs (default 1001); --min-count and
    /// --min-count-fraction act as in base coordinates. Read correction and general repeat
    /// unzipping are not implemented.
    #[arg(long, conflicts_with_all = ["multichrom_graph", "auto_k_depth"])]
    rust_hpc: bool,
    #[arg(long)]
    threads: Option<usize>,
    /// Output assembly graph (GFA, exact nM overlaps).
    #[arg(long)]
    out: PathBuf,
    #[arg(long)]
    out_json: PathBuf,
}

fn assemble(args: AssembleArgs) -> Result<()> {
    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()?;
    }
    let coordinates = if args.rust_hpc || args.multichrom_graph {
        assemble::Coordinates::HomopolymerRuns
    } else {
        assemble::Coordinates::Bases
    };
    let automatic_k = args.k.is_none() && args.s != 0;
    let report = assemble::run_engine(
        args.reads,
        assemble::EngineParams {
            cleaning_audit_only: args.cleaning_audit_only,
            coordinates,
            k: args.k,
            s: args.s,
            min_count: args.min_count,
            min_count_fraction: args.min_count_fraction,
            nuclear_share: args.nuclear_share,
            auto_k_depth: args.auto_k_depth.unwrap_or(assemble::AUTO_K_NODE_DEPTH),
            context: args.read_context.as_deref(),
        },
        &args.out,
        &args.out_json,
    )?;
    let r = match report {
        assemble::EngineReport::HomopolymerRuns(r) => {
            eprintln!(
                "[ovasm] homopolymer-run coordinates: {} unitigs, {} bp, {} links; {:.1}s (experimental)",
                r.unitigs, r.total_length, r.links, r.elapsed_seconds
            );
            return Ok(());
        }
        assemble::EngineReport::Bases(r) => r,
    };
    if automatic_k {
        let tried: Vec<String> = r
            .k_tried
            .iter()
            .map(|(k, d)| format!("k {k}: node depth {d:.0}"))
            .collect();
        eprintln!("[ovasm] automatic k {} ({})", r.k, tried.join(", "));
    }
    eprintln!(
        "[ovasm] {} reads; {} sampled k-mers, {} distinct; depth ~{:.0}, min count {}; {} nodes, {} edges -> {} unitigs ({} bp, N50 {}), {} links; {:.1}s",
        r.reads, r.sampled_kmers, r.distinct_nodes, r.node_depth, r.min_count, r.nodes_kept,
        r.edges_kept, r.unitigs, r.total_length, r.n50, r.links, r.elapsed_seconds
    );
    Ok(())
}

#[derive(clap::Args)]
struct PrimaryArgs {
    /// Optional graph-bound U1 report; constrain only fully sampled repeat ends.
    #[arg(long)]
    phases: Option<PathBuf>,
    #[arg(long)]
    gfa: PathBuf,
    #[arg(long)]
    out: PathBuf,
    #[arg(long)]
    out_json: PathBuf,
    /// Maximum segment copies inferred from relative depth.
    #[arg(long, default_value_t = 10)]
    max_copy: usize,
    /// Search expansions per connected component. Exceeding it fails without a primary FASTA.
    #[arg(long, default_value_t = 1_000_000)]
    max_steps: usize,
    /// Prefer a circular path when it covers this fraction of the longest linear path.
    #[arg(long, default_value_t = 0.9)]
    circular_fraction: f64,
}

fn primary(args: PrimaryArgs) -> Result<()> {
    let report = primary::run_phased(
        &args.gfa,
        &args.out,
        &args.out_json,
        args.max_copy,
        args.max_steps,
        args.circular_fraction,
        args.phases.as_deref(),
    )?;
    eprintln!(
        "[ovasm] primary: {} assembly records, {} bp; {} unplaced segments",
        report.records.len(),
        report.total_length,
        report.unplaced_segments
    );
    Ok(())
}

#[derive(clap::Args)]
struct UnifyArgs {
    /// Assembly graph from any backend (links must be 0M, nM or *).
    #[arg(long)]
    gfa: PathBuf,
    /// `ovasm check evidence` report for this graph: read support becomes `ev:i` on links.
    #[arg(long)]
    evidence: Option<PathBuf>,
    /// `ovasm solve linearize` report: the best path (and near-tied alternatives) become P lines.
    #[arg(long)]
    linearization: Option<PathBuf>,
    /// Sample name used for segment names and the PanSN sample field.
    #[arg(long)]
    sample: String,
    #[arg(long, default_value = "unknown")]
    backend: String,
    #[arg(long, default_value = "mitochondrion", value_parser = ["mitochondrion", "plastid"])]
    organelle: String,
    /// PanSN molecule field (default: mt for mitochondrion, pt for plastid).
    #[arg(long)]
    molecule: Option<String>,
    /// Alternatives written as extra haplotypes when the linearization is not decisive.
    #[arg(long, default_value_t = 3)]
    max_alternatives: usize,
    /// Output OV-GFA.
    #[arg(long)]
    out: PathBuf,
    /// Output JSON report.
    #[arg(long)]
    out_json: PathBuf,
}

#[derive(clap::Args)]
struct PanelArgs {
    /// A set of reads (FASTQ/FASTA, optionally .gz); repeat the flag for several sets, each is
    /// called on its own.
    #[arg(long = "reads", required = true)]
    reads: Vec<PathBuf>,
    /// Output JSON report.
    #[arg(long)]
    out_json: PathBuf,
    /// Protein panel FASTA (headers `class|gene|species`, class mito or plastid); default: the
    /// panel built into the binary.
    #[arg(long)]
    panel: Option<PathBuf>,
    /// Amino-acid k-mer length.
    #[arg(long, default_value_t = 6)]
    k: usize,
    /// Shortest stop-free stretch (residues) that is looked up.
    #[arg(long, default_value_t = 80)]
    min_orf: usize,
    /// Shared k-mers a stretch needs with its best gene.
    #[arg(long, default_value_t = 5)]
    min_hits: u32,
    /// Stretches a gene needs to count as present.
    #[arg(long, default_value_t = 2)]
    min_gene_segments: u64,
    /// Distinct genes a class needs for a call.
    #[arg(long, default_value_t = 5)]
    min_genes: usize,
    /// Gene stretches per read a class needs for a call.
    #[arg(long, default_value_t = 0.3)]
    min_density: f64,
    /// Leave panel proteins of species starting with this name out (repeatable; for testing
    /// against a genus the panel holds).
    #[arg(long = "exclude-species")]
    exclude: Vec<String>,
    #[arg(long)]
    threads: Option<usize>,
}

#[derive(clap::Args)]
struct ComponentsArgs {
    /// Complete original-read competing-context evidence (ovasm.read-origins.v1 JSON).
    #[arg(long)]
    read_origins: Option<PathBuf>,
    /// Target candidates after abstaining on components with competing origins.
    #[arg(long)]
    out_assigned_gfa: Option<PathBuf>,
    /// Preserved components whose target source cannot be assigned.
    #[arg(long)]
    out_unassigned_gfa: Option<PathBuf>,
    /// The assembly graph (GFA with sequences; depth from dp, DP, KC or RC tags).
    #[arg(long)]
    gfa: PathBuf,
    /// The organelle the graph was assembled for: mitochondrion or plastid.
    #[arg(long)]
    organelle: String,
    /// Output GFA: the kept components, lines unchanged.
    #[arg(long)]
    out: PathBuf,
    /// Output JSON report: every component with its length, depth, genes and the reason it
    /// was kept or removed.
    #[arg(long)]
    out_json: PathBuf,
    /// A component without genes is kept when its depth is within this factor of the
    /// reference component's (the one with the most target genes), either way.
    #[arg(long, default_value_t = 2.0)]
    depth_factor: f64,
    /// A component with target genes is kept unless the reference is more than this many
    /// times deeper.
    #[arg(long, default_value_t = 2.5)]
    gene_depth_factor: f64,
    /// Genes of the other organelle (and at least twice the target genes) that make a
    /// component the other organelle's.
    #[arg(long, default_value_t = 5)]
    min_other_genes: usize,
    /// Exclude an independent component with other-organelle genes and no target gene.
    /// Components lacking both classes of genes still use the depth window.
    #[arg(long)]
    reject_other_only: bool,
    /// Protein panel FASTA (headers `class|gene|species`); default: the built-in panel.
    #[arg(long)]
    panel: Option<PathBuf>,
    /// Amino-acid k-mer length.
    #[arg(long, default_value_t = 6)]
    k: usize,
    /// Shortest stop-free stretch (residues) that is looked up.
    #[arg(long, default_value_t = 80)]
    min_orf: usize,
    /// Shared k-mers a stretch needs with its best gene.
    #[arg(long, default_value_t = 5)]
    min_hits: u32,
    /// Leave panel proteins of species starting with this name out (repeatable).
    #[arg(long = "exclude-species")]
    exclude: Vec<String>,
}

fn components(args: ComponentsArgs) -> Result<()> {
    let panel = match &args.panel {
        Some(f) => panel::Panel::parse(
            &std::fs::read_to_string(f).with_context(|| format!("cannot read {}", f.display()))?,
            args.k,
            &args.exclude,
        )?,
        None => panel::Panel::embedded(args.k, &args.exclude)?,
    };
    let report = components::run(
        &args.gfa,
        &args.organelle,
        &panel,
        &components::ComponentParams {
            read_origins: args.read_origins,
            out_assigned_gfa: args.out_assigned_gfa,
            out_unassigned_gfa: args.out_unassigned_gfa,
            depth_factor: args.depth_factor,
            gene_depth_factor: args.gene_depth_factor,
            min_other_genes: args.min_other_genes,
            reject_other_only: args.reject_other_only,
            k: args.k,
            min_orf: args.min_orf,
            min_hits: args.min_hits,
        },
        &args.out,
        &args.out_json,
    )?;
    eprint!("{}", components::summary(&report));
    Ok(())
}

fn panel(args: PanelArgs) -> Result<()> {
    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()?;
    }
    let report = panel::run(
        &args.reads,
        args.panel.as_deref(),
        &panel::PanelParams {
            k: args.k,
            min_orf: args.min_orf,
            min_hits: args.min_hits,
            min_gene_segments: args.min_gene_segments,
            min_genes: args.min_genes,
            min_density: args.min_density,
        },
        &args.exclude,
        &args.out_json,
    )?;
    eprint!("{}", panel::summary(&report));
    Ok(())
}

#[derive(clap::Args)]
struct PanArgs {
    /// Sample OV-GFA (repeatable).
    #[arg(long)]
    gfa: Vec<PathBuf>,
    /// File listing sample OV-GFAs, one path per line (with or instead of --gfa).
    #[arg(long)]
    gfa_list: Option<PathBuf>,
    /// Shared k-mer length that joins samples' bases (odd, at most 63).
    #[arg(long, default_value_t = 63)]
    k: usize,
    /// Output pangenome GFA.
    #[arg(long)]
    out: PathBuf,
    /// Output JSON report.
    #[arg(long)]
    out_json: PathBuf,
}

fn pan(args: PanArgs) -> Result<()> {
    let mut inputs = args.gfa.clone();
    if let Some(list) = &args.gfa_list {
        let text = std::fs::read_to_string(list)
            .with_context(|| format!("cannot read {}", list.display()))?;
        inputs.extend(
            text.lines()
                .map(str::trim)
                .filter(|l| !l.is_empty() && !l.starts_with('#'))
                .map(PathBuf::from),
        );
    }
    if inputs.is_empty() {
        bail!("give sample graphs with --gfa or --gfa-list");
    }
    let report = pan::run(
        &inputs,
        &pan::PanParams { k: args.k },
        &args.out,
        &args.out_json,
    )?;
    eprintln!(
        "[ovasm] {} samples, {} bp -> {} pieces, {} bp ({} links); in every sample: {} pieces, {} bp; \
         verified {} segments, {} links, {} paths; {:.1}s",
        report.samples.len(),
        report.sample_bases,
        report.pieces,
        report.pan_bases,
        report.links,
        report.core_pieces,
        report.core_bases,
        report.verified_segments,
        report.verified_links,
        report.verified_paths,
        report.elapsed_seconds
    );
    Ok(())
}

fn unify(args: UnifyArgs) -> Result<()> {
    let molecule = args.molecule.clone().unwrap_or_else(|| {
        if args.organelle == "plastid" {
            "pt"
        } else {
            "mt"
        }
        .to_string()
    });
    let report = unify::run(
        &args.gfa,
        args.evidence.as_deref(),
        args.linearization.as_deref(),
        &unify::UnifyParams {
            sample: &args.sample,
            backend: &args.backend,
            organelle: &args.organelle,
            molecule: &molecule,
            max_alternatives: args.max_alternatives,
        },
        &args.out,
        &args.out_json,
    )?;
    eprintln!(
        "[ovasm] {} segments -> {} pieces ({} cut, {} bp of overlap removed); {} links; paths: {}",
        report.segments_in,
        report.pieces_out,
        report.segments_cut,
        report.overlap_bp_removed,
        report.links_out,
        if report.paths.is_empty() {
            "none".into()
        } else {
            report.paths.join(", ")
        }
    );
    Ok(())
}

#[derive(clap::Args)]
struct LinearizeArgs {
    /// Optional graph-bound U1 phase table; same conservative constraints as primary.
    #[arg(long)]
    phases: Option<PathBuf>,
    /// The assembly graph given to `ovasm check evidence`.
    #[arg(long)]
    gfa: PathBuf,
    /// JSON report written by `ovasm check evidence` for that graph.
    #[arg(long)]
    evidence: PathBuf,
    /// Output FASTA with the best linearization.
    #[arg(long)]
    out_fasta: PathBuf,
    /// Output JSON report (best path, ranked alternatives, unplaced repeats).
    #[arg(long)]
    out_json: PathBuf,
    /// Output FASTA for partial paths: what the reads join in components no cover solves, and
    /// required anchors without bridges. Not molecules (records are tagged complete=false).
    /// Default: beside --out-fasta, with `.partial` before the extension. Written only when
    /// there are partial paths; a file from an earlier run is removed otherwise.
    #[arg(long)]
    out_partial_fasta: Option<PathBuf>,
    /// Ignore phased paths and links supported by fewer reads.
    #[arg(long, default_value_t = 2)]
    min_reads: u64,
    /// Hard search budget per cover size; omitted enables one bounded adaptive retry on small graphs.
    #[arg(long)]
    max_steps: Option<u64>,
    #[arg(long, default_value_t = 10)]
    max_alternatives: usize,
    /// Anchors shallower than this fraction of the median depth are optional (low-frequency
    /// configurations the representative path may skip); without depth tags, half this fraction
    /// of the median bridge support.
    #[arg(long, default_value_t = 0.5)]
    optional_below: f64,
    /// Cover-size limit per component. 0 (automatic) uses the existing four-cover search
    /// envelope raised to the necessary bound from graph path ends. A positive value is
    /// an explicit limit; omitted cover sizes remain reported and preclude decisive output.
    #[arg(long, default_value_t = 0)]
    max_molecules: usize,
}

fn linearize(args: LinearizeArgs) -> Result<()> {
    let report = linearize::run_phased(
        &args.gfa,
        &args.evidence,
        linearize::LinearizeParams {
            min_reads: args.min_reads,
            max_steps: args.max_steps.unwrap_or(5_000_000),
            adaptive_steps: args.max_steps.is_none(),
            max_alternatives: args.max_alternatives,
            optional_below: args.optional_below,
            max_molecules: args.max_molecules,
        },
        &args.out_fasta,
        args.out_partial_fasta.as_deref(),
        &args.out_json,
        args.phases.as_deref(),
    )?;
    eprintln!(
        "[ovasm] {} anchors, {} bridges, {} components, {} candidates; {} molecules{}",
        report.anchors.len(),
        report.bridges,
        report.components.len(),
        report.candidates,
        report.molecules,
        if report.search_truncated {
            " (search truncated)"
        } else {
            ""
        }
    );
    for c in &report.components {
        match &c.best {
            Some(cover) => {
                for m in &cover.molecules {
                    eprintln!(
                        "[ovasm]   {} bp ({}): {}",
                        m.length,
                        if m.circular { "circular" } else { "linear" },
                        m.path
                    )
                }
            }
            None if c.skipped => eprintln!(
                "[ovasm]   skipped shallow component: {}",
                c.anchors.join(" ")
            ),
            None => eprintln!(
                "[ovasm]   {} anchors: no cover of {}..={} molecules{}; path ends {}",
                c.anchors.len(),
                c.min_molecules,
                c.searched_max_molecules,
                if c.search_truncated {
                    " (search truncated)"
                } else {
                    ""
                },
                c.path_ends.join(" ")
            ),
        }
    }
    if !report.unbridged_anchors.is_empty() {
        eprintln!(
            "[ovasm] warning: anchors without bridges, not placed: {}",
            report.unbridged_anchors.join(" ")
        );
    }
    if let Some(fasta) = &report.partial_fasta {
        eprintln!(
            "[ovasm] {} partial paths ({} bp; not molecules, nothing is solved by them) -> {}",
            report.partial_paths.len(),
            report.partial_paths.iter().map(|m| m.length).sum::<usize>(),
            fasta
        );
        for m in &report.partial_paths {
            eprintln!(
                "[ovasm]   partial {} bp ({}): {}",
                m.length,
                if m.circular { "circular" } else { "linear" },
                m.path
            );
        }
    }
    if !report.repeats_not_placed.is_empty() {
        eprintln!(
            "[ovasm] warning: repeats not placed: {}",
            report.repeats_not_placed.join(", ")
        );
    }
    Ok(())
}

#[derive(clap::Args)]
struct EvidenceArgs {
    /// Count exact overlap-edge words separately from long-span repeat phasing.
    #[arg(long)]
    short_read_adjacency: bool,
    /// Assembly graph with segment sequences and 0M links.
    #[arg(long)]
    gfa: PathBuf,
    /// Long reads (FASTQ/FASTA, optionally .gz); repeat for several files.
    #[arg(long = "reads", required = true)]
    reads: Vec<PathBuf>,
    /// Output JSON report.
    #[arg(long)]
    out: PathBuf,
    #[arg(long, default_value_t = 21)]
    k: usize,
    /// Junction sequence taken from each side of a link (bp).
    #[arg(long, default_value_t = 2000)]
    flank: usize,
    /// A supporting read must extend at least this far into both segments (bp).
    #[arg(long, default_value_t = 500)]
    min_anchor: usize,
    /// Chain k-mer hits whose diagonals differ by at most this (bp); raise for noisy reads.
    #[arg(long, default_value_t = 30)]
    diag_tol: i64,
    #[arg(long, default_value_t = 10)]
    min_side_hits: usize,
    /// Allowed deviation of joint spacing from the middle segment's length (fraction).
    #[arg(long, default_value_t = 0.05)]
    spacing_tol: f64,
    /// Bootstrap replicates for minority-share intervals (0 disables).
    #[arg(long, default_value_t = 1000)]
    bootstrap: usize,
    #[arg(long, default_value_t = 0)]
    seed: u64,
    /// Longest topology-only repeat inference without depth; depth-backed repeats have no cap.
    #[arg(long, default_value_t = 20_000)]
    max_repeat_len: usize,
    #[arg(long)]
    threads: Option<usize>,
}

#[derive(clap::Args)]
struct UnzipArgs {
    /// Audit-only cleaning JSON from this run; attach original-read adjacency and phase counts.
    #[arg(long)]
    cleaning_audit_json: Option<PathBuf>,
    /// Optional unchanged graph plus supported entry-repeat-exit P records.
    #[arg(long)]
    out_gfa: Option<PathBuf>,
    /// Optional graph with fully resolved one-to-one repeats split into context copies.
    #[arg(long)]
    materialize_gfa: Option<PathBuf>,
    /// Assembly graph with segment sequences (GFA).
    #[arg(long)]
    gfa: PathBuf,
    /// Original long reads (FASTQ/FASTA, optionally .gz); repeat for several files.
    #[arg(long = "reads", required = true)]
    reads: Vec<PathBuf>,
    /// Output JSON report.
    #[arg(long)]
    out_json: PathBuf,
    /// A phase counts as supported with at least this many distinct reads. Give the value the
    /// assembly resolved (its report's min_count); 3 is the engine's floor.
    #[arg(long, default_value_t = 3)]
    min_count: u32,
    /// Alignment settings of `evidence`, at its HiFi preset.
    #[arg(long, default_value_t = 21)]
    k: usize,
    #[arg(long, default_value_t = 2000)]
    flank: usize,
    #[arg(long, default_value_t = 500)]
    min_anchor: usize,
    #[arg(long, default_value_t = 30)]
    diag_tol: i64,
    #[arg(long, default_value_t = 10)]
    min_side_hits: usize,
    #[arg(long, default_value_t = 0.05)]
    spacing_tol: f64,
    #[arg(long)]
    threads: Option<usize>,
}

fn unzip(args: UnzipArgs) -> Result<()> {
    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()?;
    }
    if !(1..=kmer::MAX_K).contains(&args.k) {
        bail!("--k must be between 1 and {}", kmer::MAX_K);
    }
    let report = unzip::run_audited(
        &args.gfa,
        args.reads,
        evidence::EvidenceParams {
            short_read_adjacency: false,
            k: args.k,
            flank: args.flank,
            min_anchor: args.min_anchor,
            diag_tol: args.diag_tol,
            min_side_hits: args.min_side_hits,
            spacing_tol: args.spacing_tol,
            bootstrap: 0,
            seed: 0,
            max_repeat_len: 20_000,
        },
        args.min_count,
        &args.out_json,
        args.cleaning_audit_json.as_deref(),
    )?;
    if args.out_gfa.is_some() || args.materialize_gfa.is_some() {
        let exported = phase::export(
            &args.gfa,
            &args.out_json,
            args.out_gfa.as_deref(),
            args.materialize_gfa.as_deref(),
        )?;
        let mut json = serde_json::to_value(&report)?;
        json["export"] = serde_json::to_value(exported)?;
        std::fs::write(&args.out_json, serde_json::to_string_pretty(&json)?)?;
    }
    let count = |class: &str| {
        report
            .repeat_nodes
            .iter()
            .filter(|n| n.classification == class)
            .count()
    };
    eprintln!(
        "[ovasm] unzip: {} repeat nodes (resolved {}, partial {}, mixed {}, conflicted {}, unresolved {}); {} of {} reads complete a traversal; {:.1}s",
        report.repeat_nodes.len(),
        count("resolved"),
        count("partial"),
        count("mixed"),
        count("conflicted"),
        count("unresolved"),
        report.reads_with_complete_traversal,
        report.unique_records,
        report.elapsed_seconds
    );
    Ok(())
}

fn evidence(args: EvidenceArgs) -> Result<()> {
    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()?;
    }
    if !(1..=kmer::MAX_K).contains(&args.k) {
        bail!("--k must be between 1 and {}", kmer::MAX_K);
    }
    let report = evidence::run(
        &args.gfa,
        args.reads,
        evidence::EvidenceParams {
            short_read_adjacency: args.short_read_adjacency,
            k: args.k,
            flank: args.flank,
            min_anchor: args.min_anchor,
            diag_tol: args.diag_tol,
            min_side_hits: args.min_side_hits,
            spacing_tol: args.spacing_tol,
            bootstrap: args.bootstrap,
            seed: args.seed,
            max_repeat_len: args.max_repeat_len,
        },
        &args.out,
    )?;
    eprintln!(
        "[ovasm] links supported {}/{}; branches {}; repeat pairings {}; phased paths {}; {} of {} reads informative; {:.1}s",
        report.links_supported, report.links.len(), report.branches.len(), report.repeats.len(),
        report.phased_paths.len(),
        report.reads_with_support, report.reads_scanned, report.elapsed_seconds
    );
    Ok(())
}

#[derive(clap::Args)]
struct DiscoverArgs {
    /// Read files (FASTQ/FASTA, optionally .gz).
    #[arg(long = "reads", required = true)]
    reads: Vec<PathBuf>,
    #[arg(long)]
    out: PathBuf,
    #[arg(long, default_value_t = 21)]
    k: usize,
    /// Keep 1/scale of k-mers in the depth sketch. Long reads: 64; short reads need ~2-4.
    #[arg(long, default_value_t = 64)]
    scale: u64,
    /// Gigabases read to build the depth sketch; 0 (default) reads the whole input, which a
    /// prefix only stands for when the reads come in random order.
    #[arg(long, default_value_t = 0.0)]
    sketch_gb: f64,
    /// Organelle threshold as a multiple of the nuclear depth peak.
    #[arg(long, default_value_t = 3.0)]
    depth_fold: f64,
    #[arg(long, default_value_t = 20)]
    min_sampled: usize,
    #[arg(long, default_value_t = 0.9)]
    min_high_frac: f64,
    /// Maximum (Q3-Q1)/median of sampled depths along a read (uniformity; IR/SC or repeat
    /// boundaries legitimately give 1x/2x mixtures, so keep this loose).
    #[arg(long, default_value_t = 1.5)]
    max_spread: f64,
    /// Minimum distinct/total sampled k-mers in a read (rejects tandem arrays such as rDNA).
    #[arg(long, default_value_t = 0.8)]
    min_distinct: f64,
    #[arg(long, default_value_t = 20)]
    min_cluster_reads: usize,
    #[arg(long)]
    target_depth: Option<f64>,
    #[arg(long, default_value_t = 0)]
    salt: u64,
    #[arg(long)]
    threads: Option<usize>,
    /// Write per-read features to profiles.tsv.
    #[arg(long)]
    dump_profiles: bool,
}

fn discover(args: DiscoverArgs) -> Result<()> {
    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()?;
    }
    if !(1..=kmer::MAX_K).contains(&args.k) {
        bail!("--k must be between 1 and {}", kmer::MAX_K);
    }
    let report = seedfree::run(
        args.reads,
        seedfree::DiscoverParams {
            k: args.k,
            scale: args.scale,
            sketch_bases: (args.sketch_gb * 1e9) as u64,
            depth_fold: args.depth_fold,
            min_sampled: args.min_sampled,
            min_high_frac: args.min_high_frac,
            max_spread: args.max_spread,
            min_distinct: args.min_distinct,
            min_cluster_reads: args.min_cluster_reads,
            target_depth: args.target_depth,
            salt: args.salt,
            dump_profiles: args.dump_profiles,
        },
        args.out,
    )?;
    for c in &report.clusters {
        eprintln!(
            "[ovasm] {}: ~{:.0}x ({} reads, {:.1} Mb); kept {} -> {}",
            c.name,
            c.estimated_depth,
            c.reads,
            c.bases as f64 / 1e6,
            c.kept_reads,
            c.output
        );
    }
    eprintln!(
        "[ovasm] organelle-like reads: {} of {}; elapsed {:.1}s; peak RSS {}",
        report.organelle_like_reads,
        report.reads_scanned,
        report.elapsed_seconds,
        report
            .peak_rss_mb
            .map_or("n/a".into(), |m| format!("{m:.0} MB"))
    );
    Ok(())
}

#[derive(clap::Args)]
struct RecruitArgs {
    /// Read files (FASTQ/FASTA, optionally .gz); repeat for several libraries.
    #[arg(long = "reads", required = true)]
    reads: Vec<PathBuf>,
    /// Extend a target only from reads whose new k-mers are at its depth in a count of the
    /// whole library (4 GiB, after the first pass): keeps short-read extension out of the
    /// nucleus. `run` uses it for short reads.
    #[arg(long)]
    kmer_gate: bool,
    /// Downsample by the first word of each read name, so a pair's mates stay together.
    #[arg(long)]
    mates_together: bool,
    /// Seed as NAME=PATH (e.g. mito=ref_mt.fa); repeat for more targets or more seed files.
    #[arg(long = "seed")]
    seeds: Vec<String>,
    /// Seed from reads as NAME=PATH: the k-mers solid in them (seen in --extend-min-count reads
    /// and within a tenth of the depth peak), e.g. corrected long reads of the same sample to
    /// recruit its short reads. Combines with --seed; one of the two is required.
    #[arg(long = "seed-reads")]
    seed_reads: Vec<String>,
    /// Output directory (per-target FASTQ, numt_candidates.fastq, recruit.json).
    #[arg(long)]
    out: PathBuf,
    /// Read platform preset: sets k, block_gap, min_hit_frac and numt_min_block defaults.
    #[arg(long, value_enum, default_value_t = Preset::Hifi)]
    preset: Preset,
    /// k-mer size (<= 31); smaller tolerates more divergence/errors. Overrides the preset.
    #[arg(long)]
    k: Option<usize>,
    #[arg(long, default_value_t = 0.8)]
    min_cover: f64,
    /// Overrides the preset.
    #[arg(long)]
    min_hit_frac: Option<f64>,
    /// Minimum organelle-like block (bp) for a partially matching read to be a NUMT candidate.
    #[arg(long)]
    numt_min_block: Option<usize>,
    /// Largest gap (bp) between bait hits still bridged into one block. Overrides the preset.
    #[arg(long)]
    block_gap: Option<usize>,
    /// A read the bait covers but with fewer of its own k-mers distinct than this share is a
    /// tandem array and is not recruited. Default: 0.5 for hifi (see recruit::HIFI_MIN_DISTINCT),
    /// 0 (off) for sr and ont.
    #[arg(long)]
    min_distinct: Option<f64>,
    /// Most scanning passes; passes after the first add k-mers from recruited reads (iterative
    /// baiting). Iteration stops as soon as no target grows by --saturation, so a close seed
    /// costs one extra pass while a distant one keeps extending until complete. Default: 20 for
    /// hifi and sr; 1 for ont, whose reads (ONT, CLR) are too noisy to extend the bait: their
    /// error k-mers are short random words that also occur in the nuclear genome, and extension
    /// runs away into nuclear reads. Correct noisy reads first, then extend from them.
    #[arg(long)]
    iterations: Option<usize>,
    /// A k-mer must occur in at least this many recruited reads to join the bait set.
    #[arg(long, default_value_t = 3)]
    extend_min_count: u32,
    /// Internal automatic HiFi pipeline policy; never override an explicit support floor.
    #[arg(long, hide = true, conflicts_with = "extend_min_count")]
    auto_shallow_edge_support: bool,
    /// New reads of a target are subsampled to this depth before bait extension.
    #[arg(long, default_value_t = 200.0)]
    extend_depth: f64,
    /// Stop extending a target whose recruited reads grew by less than this fraction (with
    /// --extend-gate, unless they grew by at least half as much as in the pass before: a large
    /// genome grows by about as many reads every pass, a smaller share each time).
    #[arg(long, default_value_t = 0.05)]
    saturation: f64,
    /// Downsample each target to about this depth (recruited bases / seed length).
    #[arg(long)]
    target_depth: Option<f64>,
    /// Downsample each target to about this read depth: the median depth of its reads in the
    /// depth sketch, so it needs no genome size (a multi-species seed library has none, and
    /// --target-depth divides by its summed length). Hifi preset only; builds the sketch.
    #[arg(long)]
    max_read_depth: Option<f64>,
    /// Downsample each target to at most about this many bases.
    #[arg(long)]
    max_bases: Option<u64>,
    /// New scratch directory: classify original inputs once, replay sequences, then restore original qualities in one final input pass.
    #[arg(long)]
    recruit_cache: Option<PathBuf>,
    /// Salt for deterministic downsampling; change it to draw an independent replicate.
    #[arg(long, default_value_t = 0)]
    salt: u64,
    /// Worker threads (default: all cores).
    #[arg(long)]
    threads: Option<usize>,
    /// Bootstrap a target whose seed is too distant to recruit its reads: once every other
    /// target's bait has stopped growing, partially matching reads that are deep along their whole length
    /// (the `find discover` test) seed it. Seed every organelle of the sample (e.g. mitochondrion and
    /// plastid) so that each target's reads are claimed by their own. Long reads (hifi) only.
    /// Experimental: when another target's seed recruits nothing, that target counts as
    /// complete and the bootstrap can run away into nuclear repeats (v2 holdout, Jasione).
    #[arg(long)]
    bootstrap: bool,
    /// Extend each target's bait only from reads that are organelle-deep along their whole
    /// length (the `find discover` test) and not shallower than a quarter of the target's own reads,
    /// so a stray seed sequence cannot lead extension into the nucleus (nor a plastid into the
    /// mitochondrion); a plastid or mitochondrion whose reads come to stand for more sequence
    /// than any such genome (1 and 20 Mb) stops extending. Use with `--saturation 0` to extend
    /// until no bait grows. Long reads (hifi) only.
    #[arg(long)]
    extend_gate: bool,
    /// The depth gate's organelle threshold, in multiples of the nuclear peak depth (default 3).
    /// Extension already lowers it for a target shallower than that (see `--extend-gate`).
    #[arg(long)]
    gate_fold: Option<f64>,
    /// The depth gate's nuclear peak (sketch depth, about the single-copy read depth), given
    /// instead of judged from the reads' median depths. Those stand two copies deep in a
    /// repeat-rich or polyploid genome (Lespedeza davurica: 56x for 27x).
    #[arg(long, conflicts_with = "genome_size")]
    nuclear_depth: Option<u32>,
    /// The nuclear genome's size, for the depth gate's nuclear peak (the input's bases over it):
    /// the length of a collapsed assembly, all subgenomes of an allopolyploid (1C) but one set of
    /// an autopolyploid; bases, or with k, m or g (1.26g).
    #[arg(long, value_parser = parse_genome_size)]
    genome_size: Option<u64>,
}

/// Bases, plain or with a k, m or g suffix (1.26g).
fn parse_bases(s: &str) -> std::result::Result<u64, String> {
    let t = s.trim().to_ascii_lowercase();
    let (num, scale) = match t.strip_suffix(['k', 'm', 'g']) {
        Some(n) => (n, match t.chars().last() { Some('k') => 1e3, Some('m') => 1e6, _ => 1e9 }),
        None => (t.as_str(), 1.0),
    };
    match num.parse::<f64>() {
        Ok(v) if v.is_finite() && v * scale >= 1.0 => Ok((v * scale).round() as u64),
        _ => Err(format!("{s}: not a size in bases (e.g. 1260000000 or 1.26g)")),
    }
}

/// A nuclear genome's size (`parse_bases`) within land plants' range, so that a forgotten
/// suffix (1260 for 1260m) or a stray exponent does not set the gate's nuclear depth to
/// thousands or to 1.
fn parse_genome_size(s: &str) -> std::result::Result<u64, String> {
    let b = parse_bases(s)?;
    if (MIN_GENOME_SIZE..=MAX_GENOME_SIZE).contains(&b) {
        Ok(b)
    } else {
        Err(format!(
            "{s}: {b} bases is no nuclear genome (land plants have about 60 Mb to 150 Gb); give the size with k, m or g (1.26g)"
        ))
    }
}

/// Bounds of `parse_genome_size`: below the smallest land-plant nuclear genomes (~60 Mb) and
/// above the largest (Paris japonica, ~150 Gb).
const MIN_GENOME_SIZE: u64 = 10_000_000;
const MAX_GENOME_SIZE: u64 = 300_000_000_000;

#[derive(clap::Args)]
struct LabelArgs {
    /// Read clusters (e.g. `ovasm find discover`'s clusterN.fastq), each assembled on its own;
    /// repeat for every cluster.
    #[arg(long = "reads", required = true)]
    reads: Vec<PathBuf>,
    /// Organelle protein database: FASTA with headers organelle|gene|source.
    #[arg(long)]
    genes: PathBuf,
    /// Output directory (per-cluster graphs, seed_<organelle>.fasta, labels.json).
    #[arg(long)]
    out: PathBuf,
    #[arg(long)]
    threads: Option<usize>,
}

fn label(args: LabelArgs) -> Result<()> {
    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()?;
    }
    let r = label::run(
        &args.reads,
        &args.genes,
        label::LabelParams::default(),
        &args.out,
    )?;
    for (o, name) in label::ORGANELLES.iter().enumerate() {
        eprintln!(
            "[ovasm] {name}: {} bp of seed from {} clusters",
            r.seed_bp[o],
            r.clusters.iter().filter(|c| c.seed_bp[o] > 0).count()
        );
    }
    Ok(())
}

#[derive(Clone, Copy, clap::ValueEnum)]
enum Preset {
    /// PacBio HiFi (>Q20): exact 21-mers, short gap bridging.
    Hifi,
    /// Raw Oxford Nanopore / CLR: shorter k-mers, wide gap bridging, low hit-fraction floor.
    Ont,
    /// Illumina short reads: bridge SNV-sized gaps; NUMT flagging needs long reads.
    Sr,
}

impl Preset {
    /// (k, block_gap, min_hit_frac, numt_min_block)
    fn defaults(self) -> (usize, usize, f64, usize) {
        match self {
            Preset::Hifi => (21, 50, 0.2, 1000),
            Preset::Ont => (15, 200, 0.05, 2000),
            Preset::Sr => (21, 40, 0.2, usize::MAX),
        }
    }
}

fn recruit(args: RecruitArgs) -> Result<()> {
    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()?;
    }
    let (pk, pgap, phit, pnumt) = args.preset.defaults();
    let k = args.k.unwrap_or(pk);
    if !(1..=kmer::MAX_K).contains(&k) {
        bail!("--k must be between 1 and {}", kmer::MAX_K);
    }
    if args.bootstrap && !matches!(args.preset, Preset::Hifi) {
        bail!("--bootstrap needs the hifi preset: its depth gate judges whole long reads");
    }
    if args.extend_gate && !matches!(args.preset, Preset::Hifi) {
        bail!("--extend-gate needs the hifi preset: its depth gate judges whole long reads");
    }
    if args.gate_fold.is_some_and(|f| f.is_nan() || f <= 0.0) {
        bail!("--gate-fold must be positive");
    }
    if args.gate_fold.is_some()
        && !(args.extend_gate || args.bootstrap || args.max_read_depth.is_some())
    {
        bail!(
            "--gate-fold sets the depth gate's organelle threshold: give --extend-gate, --bootstrap or --max-read-depth"
        );
    }
    if args.nuclear_depth == Some(0) {
        bail!("--nuclear-depth must be positive");
    }
    if (args.nuclear_depth.is_some() || args.genome_size.is_some())
        && !(args.extend_gate || args.bootstrap || args.max_read_depth.is_some())
    {
        bail!(
            "--nuclear-depth and --genome-size set the depth gate's nuclear peak: give \
             --extend-gate, --bootstrap or --max-read-depth"
        );
    }
    if args.max_read_depth.is_some() && !matches!(args.preset, Preset::Hifi) {
        bail!("--max-read-depth needs the hifi preset: its depth sketch judges whole long reads");
    }
    if args.max_read_depth.is_some_and(|d| d.is_nan() || d <= 0.0) || args.max_bases == Some(0)
    {
        bail!("--max-read-depth and --max-bases must be positive");
    }
    let mut bait = bait::BaitSet::new(k);
    if args.seeds.is_empty() && args.seed_reads.is_empty() {
        bail!("give --seed NAME=PATH or --seed-reads NAME=PATH");
    }
    for spec in &args.seeds {
        let (name, path) = spec
            .split_once('=')
            .with_context(|| format!("--seed must be NAME=PATH, got {spec:?}"))?;
        let t = bait.target(name)?;
        bait.add_seed_file(t, &PathBuf::from(path))?;
    }
    for spec in &args.seed_reads {
        let (name, path) = spec
            .split_once('=')
            .with_context(|| format!("--seed-reads must be NAME=PATH, got {spec:?}"))?;
        let t = bait.target(name)?;
        let n = bait.add_seed_reads(t, &PathBuf::from(path), args.extend_min_count)?;
        eprintln!("[ovasm] seed {name}: {n} solid k-mers from reads {path}");
    }
    eprintln!(
        "[ovasm] bait: {} k-mers, k={}, targets={:?}",
        bait.len(),
        bait.k,
        bait.targets.iter().map(|t| &t.name).collect::<Vec<_>>()
    );
    let report = recruit::run(
        bait,
        recruit::RecruitConfig {
            bait_growth_guard: false,
            replay_dir: args.recruit_cache,
            inputs: args.reads,
            params: recruit::ClassifyParams {
                min_distinct: args.min_distinct.unwrap_or(match args.preset {
                    Preset::Hifi => recruit::HIFI_MIN_DISTINCT,
                    Preset::Ont | Preset::Sr => 0.0,
                }),
                min_cover: args.min_cover,
                min_hit_frac: args.min_hit_frac.unwrap_or(phit),
                numt_min_block: args.numt_min_block.unwrap_or(pnumt),
                block_gap: args.block_gap.unwrap_or(pgap),
            },
            iterations: args.iterations.unwrap_or(match args.preset {
                Preset::Ont => 1,
                Preset::Hifi | Preset::Sr => 20,
            }),
            extend_min_count: args.extend_min_count,
            auto_shallow_edge_support: args.auto_shallow_edge_support,
            extend_depth: args.extend_depth,
            saturation: args.saturation,
            target_depth: args.target_depth,
            salt: args.salt,
            out_dir: args.out,
            bootstrap: args.bootstrap,
            extend_gate: args.extend_gate,
            max_read_depth: args.max_read_depth,
            max_bases: args.max_bases,
            mates_together: args.mates_together,
            kmer_gate: args.kmer_gate,
            gate_fold: args.gate_fold,
            nuclear_depth: args.nuclear_depth,
            genome_size: args.genome_size,
        },
    )?;
    for o in &report.outputs {
        let depth = match (o.read_depth, o.implied_size) {
            (Some(d), Some(s)) => format!(" (read depth {d:.0}x, standing for {:.0} kb)", s as f64 / 1e3),
            _ => String::new(),
        };
        eprintln!(
            "[ovasm] {}: recruited {} reads ({:.1} Mb, ~{:.0}x){depth}; {}; kept {} reads (fraction {:.3}) -> {}",
            o.name, o.recruited_reads, o.recruited_bases as f64 / 1e6, o.estimated_depth,
            o.stop_reason, o.kept_reads, o.sampling_fraction, o.output
        );
    }
    eprintln!(
        "[ovasm] NUMT candidates: {}; elapsed {:.1}s; peak RSS {}",
        report.numt_candidates,
        report.elapsed_seconds,
        report
            .peak_rss_mb
            .map_or("n/a".into(), |m| format!("{m:.0} MB"))
    );
    Ok(())
}

/// Run one subcommand; `pipeline` runs its stages through this too.
fn dispatch(command: Command) -> Result<()> {
    match command {
        Command::Run(RunCommand::Short(mut args)) => {
            args.read_type = pipeline::ReadType::Sr;
            run(args)
        }
        Command::Run(RunCommand::Hifi(mut args)) => {
            args.read_type = pipeline::ReadType::Hifi;
            run(args)
        }
        Command::Run(RunCommand::Ont(args)) => longread::run_ont(args),
        Command::Run(RunCommand::Hybrid(args)) => longread::run_hybrid(args),
        Command::Find(FindCommand::Recruit(args)) => recruit(args),
        Command::Find(FindCommand::Discover(args)) => discover(args),
        Command::Find(FindCommand::Classify(args)) => panel(args),
        Command::Find(FindCommand::Label(args)) => label(args),
        Command::Assemble(AssembleCommand::Correct(args)) => correct(args),
        Command::Assemble(AssembleCommand::Graph(args)) => assemble(args),
        Command::Assemble(AssembleCommand::Filter(args)) => components(args),
        Command::Assemble(AssembleCommand::Unify(args)) => unify(args),
        Command::Check(CheckCommand::Evidence(args)) => evidence(args),
        Command::Check(CheckCommand::Unzip(args)) => unzip(args),
        Command::Check(CheckCommand::Assess(args)) => assess(args),
        Command::Solve(SolveCommand::Linearize(args)) => linearize(args),
        Command::Solve(SolveCommand::Primary(args)) => primary(args),
        Command::Solve(SolveCommand::Chromosomes(args)) => chromosomes(args),
        Command::Pan(PanCommand::Import(args)) => {
            pan_external::run(&pan_external::Config { input: args.input, extra_paths: args.extra_paths, out: args.out, sample: args.sample })?;
            Ok(())
        }
        Command::Pan(PanCommand::Build(args)) => build_external_pan(args),
        Command::Pan(PanCommand::Merge(args)) => pan(args),
        Command::Pan(PanCommand::Reconstruct(args)) => reconstruct(args),
        Command::Pan(PanCommand::Identify(args)) => identify(args),
        Command::Pan(PanCommand::Sharing(args)) => {
            if !(1..=kmer::MAX_K).contains(&args.k) {
                bail!("--k must be between 1 and {}", kmer::MAX_K);
            }
            retained::cross_organelle(&args.assembly, &args.out, args.k)?;
            Ok(())
        }
    }
}

fn main() -> Result<()> {
    dispatch(Cli::parse().command)
}

#[cfg(test)]
mod tests {
    use super::{parse_bases, parse_genome_size};
    #[test]
    fn explicit_support_floor_conflicts_with_automatic_shallow_policy() {
        use clap::Parser;
        let base = ["ovasm", "find", "recruit", "--reads", "reads.fa", "--out", "out"];
        for floor in ["3", "10"] {
            let mut argv = base.to_vec();
            argv.extend(["--extend-min-count", floor]);
            let parsed = super::Cli::try_parse_from(&argv).unwrap();
            if let super::Command::Find(super::FindCommand::Recruit(args)) = parsed.command {
                assert!(!args.auto_shallow_edge_support);
                assert_eq!(args.extend_min_count, floor.parse::<u32>().unwrap());
            } else { panic!("wrong command"); }
            argv.push("--auto-shallow-edge-support");
            assert!(super::Cli::try_parse_from(argv).is_err());
        }
        let mut automatic = base.to_vec();
        automatic.push("--auto-shallow-edge-support");
        let parsed = super::Cli::try_parse_from(automatic).unwrap();
        if let super::Command::Find(super::FindCommand::Recruit(args)) = parsed.command {
            assert!(args.auto_shallow_edge_support);
            assert_eq!(args.extend_min_count, 3);
        } else { panic!("wrong command"); }
    }


    #[test]
    fn genome_sizes_parse_with_or_without_a_suffix() {
        assert_eq!(parse_bases("1260000000"), Ok(1_260_000_000));
        assert_eq!(parse_bases("1.26g"), Ok(1_260_000_000));
        assert_eq!(parse_bases("1.26G"), Ok(1_260_000_000));
        assert_eq!(parse_bases("450m"), Ok(450_000_000));
        assert_eq!(parse_bases("150k"), Ok(150_000));
        assert!(parse_bases("0").is_err());
        assert!(parse_bases("-3g").is_err());
        assert!(parse_bases("big").is_err());
    }

    #[test]
    fn genome_sizes_outside_land_plants_are_refused() {
        assert_eq!(parse_genome_size("1.26g"), Ok(1_260_000_000));
        assert!(parse_genome_size("1260").is_err(), "a forgotten suffix");
        assert!(parse_genome_size("150k").is_err());
        assert!(parse_genome_size("1e30g").is_err());
    }
}
