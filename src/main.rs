//! ovasm: OrganelleVerse organelle assembler (run, recruit, discover, assemble, components, evidence,
//! linearize, unify, pan, identify).

mod assemble;
mod bait;
mod components;
mod correct;
mod evidence;
mod identify;
mod kmer;
mod label;
mod linearize;
mod pan;
mod panel;
mod pipeline;
mod recruit;
mod replicate;
mod seeddb;
mod seedfree;
mod unify;

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
    /// Reads in, organelle graph and molecules out: recruit, assemble (automatic k, smaller k
    /// until the molecules explain half of the graph), evidence, linearize, unify; --mode
    /// standard adds downsampled replicates. The other subcommands are its stages.
    Run(RunArgs),
    /// Recruit mitochondrial/plastid reads from whole-genome sequencing data in one pass.
    Recruit(RecruitArgs),
    /// Seed-free discovery: find organelle reads by k-mer depth, split them into depth clusters.
    Discover(DiscoverArgs),
    /// Read support for every link, branch and repeat pairing of an assembly graph (GFA).
    Evidence(EvidenceArgs),
    /// Representative linear (or circular) sequence through read-supported phased paths.
    Linearize(LinearizeArgs),
    /// Rewrite a backend's assembly graph as OV-GFA v1: blunt links, PanSN paths, common tags.
    Unify(UnifyArgs),
    /// Merge samples' OV-GFA graphs into one pangenome graph, keeping every configuration path.
    Pan(PanArgs),
    /// Say which organelle each set of reads (e.g. a `discover` cluster) comes from, by
    /// conserved mitochondrial and plastid proteins; no reference genome needed.
    Panel(PanelArgs),
    /// Assemble HiFi organelle reads on a closed-syncmer sparse de Bruijn graph (GFA out).
    Assemble(AssembleArgs),
    /// Self-correct noisy organelle reads (ONT, CLR) against their own solid k-mers.
    Correct(CorrectArgs),
    /// Closest reference genomes of an assembly by k-mer containment (organelle and species).
    Identify(IdentifyArgs),
    /// Assemble each read cluster and label unitigs mitochondrion or plastid by gene content;
    /// write the labelled unitigs as the sample's own seeds.
    Label(LabelArgs),
    /// Keep the connected components of an assembly graph that are the target organelle (by
    /// conserved proteins and depth); for graphs of reads found by `discover`.
    Components(ComponentsArgs),
}

#[derive(clap::Args)]
struct IdentifyArgs {
    /// Assembled sequences (FASTA), e.g. `molecules.fasta` from `ovasm linearize`.
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
    /// hifi: sparse graph at the automatic k, retried at the next smaller k of 1001, 701, 501,
    /// 351, 251, 201, 151, 101 until the molecules add up to at least half of the graph's
    /// length (when no k does, the attempt that comes closest is kept and summary.json says
    /// k_accepted: false). sr (Illumina): dense graph, k from the read length, no retry.
    #[arg(long, value_enum, default_value_t = pipeline::ReadType::Hifi)]
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
    /// components that are not the target organelle (`ovasm components`; assembly.filtered.gfa,
    /// components.json), and the later stages work on the kept ones.
    #[arg(long, value_enum, default_value_t = pipeline::RecruitBy::Seeds)]
    recruit: pipeline::RecruitBy,
    /// Downsample seed-recruited reads to about this depth (recruited bases / seed length).
    /// Right for a single close reference (the desktop runner uses 150 for a custom seed);
    /// leave unset for a multi-species seed database, whose summed length is no genome size.
    #[arg(long)]
    recruit_depth: Option<f64>,
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
    /// built-in panel. See `ovasm components`.
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
                builtin_seeds.iter().cloned().collect::<Vec<_>>().join(" and "),
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
    /// Organelle reads (e.g. from `ovasm recruit`), FASTQ/FASTA, optionally .gz.
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
    /// Organelle reads (e.g. from `ovasm recruit`), FASTQ/FASTA, optionally .gz.
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
    let r = match args.k {
        None if args.s != 0 => {
            let r = assemble::run_auto_k(
                args.reads,
                args.s,
                args.min_count,
                args.auto_k_depth.unwrap_or(assemble::AUTO_K_NODE_DEPTH),
                &args.out,
                &args.out_json,
            )?;
            let tried: Vec<String> = r
                .k_tried
                .iter()
                .map(|(k, d)| format!("k {k}: node depth {d:.0}"))
                .collect();
            eprintln!("[ovasm] automatic k {} ({})", r.k, tried.join(", "));
            r
        }
        k => {
            let k = match k {
                Some(k) => k,
                None => {
                    let (k, median) = assemble::dense_k(&args.reads)?;
                    eprintln!("[ovasm] dense mode: median read length {median} bp, k {k}");
                    k
                }
            };
            assemble::run(
                args.reads,
                assemble::AssembleParams {
                    k,
                    s: args.s,
                    min_count: args.min_count,
                },
                &args.out,
                &args.out_json,
            )?
        }
    };
    eprintln!(
        "[ovasm] {} reads; {} sampled k-mers, {} distinct; depth ~{:.0}, min count {}; {} nodes, {} edges -> {} unitigs ({} bp, N50 {}), {} links; {:.1}s",
        r.reads, r.sampled_kmers, r.distinct_nodes, r.node_depth, r.min_count, r.nodes_kept,
        r.edges_kept, r.unitigs, r.total_length, r.n50, r.links, r.elapsed_seconds
    );
    Ok(())
}

#[derive(clap::Args)]
struct UnifyArgs {
    /// Assembly graph from any backend (links must be 0M, nM or *).
    #[arg(long)]
    gfa: PathBuf,
    /// `ovasm evidence` report for this graph: read support becomes `ev:i` on links.
    #[arg(long)]
    evidence: Option<PathBuf>,
    /// `ovasm linearize` report: the best path (and near-tied alternatives) become P lines.
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
            depth_factor: args.depth_factor,
            gene_depth_factor: args.gene_depth_factor,
            min_other_genes: args.min_other_genes,
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
    /// The assembly graph given to `ovasm evidence`.
    #[arg(long)]
    gfa: PathBuf,
    /// JSON report written by `ovasm evidence` for that graph.
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
    #[arg(long, default_value_t = 5_000_000)]
    max_steps: u64,
    #[arg(long, default_value_t = 10)]
    max_alternatives: usize,
    /// Anchors shallower than this fraction of the median depth are optional (low-frequency
    /// configurations the representative path may skip); without depth tags, half this fraction
    /// of the median bridge support.
    #[arg(long, default_value_t = 0.5)]
    optional_below: f64,
    /// Most molecules (sub-circles, chromosomes sharing repeats) one connected part of the graph
    /// may be unringed into before it is reported unsolved.
    #[arg(long, default_value_t = 4)]
    max_molecules: usize,
}

fn linearize(args: LinearizeArgs) -> Result<()> {
    let report = linearize::run(
        &args.gfa,
        &args.evidence,
        linearize::LinearizeParams {
            min_reads: args.min_reads,
            max_steps: args.max_steps,
            max_alternatives: args.max_alternatives,
            optional_below: args.optional_below,
            max_molecules: args.max_molecules,
        },
        &args.out_fasta,
        args.out_partial_fasta.as_deref(),
        &args.out_json,
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
                args.max_molecules,
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
    /// New reads of a target are subsampled to this depth before bait extension.
    #[arg(long, default_value_t = 200.0)]
    extend_depth: f64,
    /// Stop extending a target whose recruited reads grew by less than this fraction.
    #[arg(long, default_value_t = 0.05)]
    saturation: f64,
    /// Downsample each target to about this depth (recruited bases / seed length).
    #[arg(long)]
    target_depth: Option<f64>,
    /// Salt for deterministic downsampling; change it to draw an independent replicate.
    #[arg(long, default_value_t = 0)]
    salt: u64,
    /// Worker threads (default: all cores).
    #[arg(long)]
    threads: Option<usize>,
    /// Bootstrap a target whose seed is too distant to recruit its reads: once every other
    /// target's bait has stopped growing, partially matching reads that are deep along their whole length
    /// (the `discover` test) seed it. Seed every organelle of the sample (e.g. mitochondrion and
    /// plastid) so that each target's reads are claimed by their own. Long reads (hifi) only.
    /// Experimental: when another target's seed recruits nothing, that target counts as
    /// complete and the bootstrap can run away into nuclear repeats (v2 holdout, Jasione).
    #[arg(long)]
    bootstrap: bool,
    /// Extend each target's bait only from reads that are organelle-deep along their whole
    /// length (the `discover` test), so a stray seed sequence cannot lead extension into the
    /// nucleus; use with `--saturation 0` to extend until no bait grows. Long reads (hifi) only.
    #[arg(long)]
    extend_gate: bool,
}

#[derive(clap::Args)]
struct LabelArgs {
    /// Read clusters (e.g. `ovasm discover`'s clusterN.fastq), each assembled on its own;
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
            inputs: args.reads,
            params: recruit::ClassifyParams {
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
            extend_depth: args.extend_depth,
            saturation: args.saturation,
            target_depth: args.target_depth,
            salt: args.salt,
            out_dir: args.out,
            bootstrap: args.bootstrap,
            extend_gate: args.extend_gate,
        },
    )?;
    for o in &report.outputs {
        eprintln!(
            "[ovasm] {}: recruited {} reads ({:.1} Mb, ~{:.0}x); kept {} reads (fraction {:.3}) -> {}",
            o.name, o.recruited_reads, o.recruited_bases as f64 / 1e6, o.estimated_depth,
            o.kept_reads, o.sampling_fraction, o.output
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
        Command::Run(args) => run(args),
        Command::Recruit(args) => recruit(args),
        Command::Discover(args) => discover(args),
        Command::Evidence(args) => evidence(args),
        Command::Linearize(args) => linearize(args),
        Command::Unify(args) => unify(args),
        Command::Assemble(args) => assemble(args),
        Command::Correct(args) => correct(args),
        Command::Pan(args) => pan(args),
        Command::Identify(args) => identify(args),
        Command::Panel(args) => panel(args),
        Command::Label(args) => label(args),
        Command::Components(args) => components(args),
    }
}

fn main() -> Result<()> {
    dispatch(Cli::parse().command)
}
