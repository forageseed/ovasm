//! Optional external pan builders. No shell, no inferred sample paths.
use std::{collections::{BTreeMap, BTreeSet}, fs, io::{Read, Write}, path::{Path, PathBuf}, process::{Command, Stdio}};
use anyhow::{bail, Context, Result};
use serde_json::{json, Value};

#[derive(Clone, Copy, Debug)]
pub enum Backend { Pggb, Minigraph }
impl Backend { fn name(self) -> &'static str { match self { Self::Pggb => "pggb", Self::Minigraph => "minigraph" } } }
pub struct Config {
    pub backend: Backend,
    /// One haploid assembly per labeled sample; minigraph uses the first as backbone.
    pub inputs: Vec<(String, PathBuf)>,
    pub out: PathBuf,
    pub threads: usize,
    pub executable: Option<PathBuf>,
    pub sample: String,
}
fn safe(s: &str) -> bool { !s.is_empty() && s.bytes().all(|b| b.is_ascii_alphanumeric() || b"_.-".contains(&b)) }
fn prepare(path: &Path, sample: &str, output: &Path, mapping: &mut String, provenance: &mut BTreeMap<String, Vec<String>>) -> Result<()> {
    let file = fs::File::open(path).with_context(|| format!("read FASTA {}", path.display()))?;
    let mut reader: Box<dyn Read> = if path.extension().is_some_and(|e| e == "gz") {
        Box::new(flate2::read::MultiGzDecoder::new(file))
    } else { Box::new(file) };
    let mut text = String::new(); reader.read_to_string(&mut text)?;
    let mut records: Vec<(String,String)>=Vec::new();
    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        if let Some(header)=line.strip_prefix('>') { records.push((header.to_string(),String::new())); }
        else {
            let record=records.last_mut().context("expected assembled FASTA, not reads")?;
            let seq=line.trim().to_ascii_uppercase();
            if !seq.bytes().all(|b| b"ACGTRYSWKMBDHVN".contains(&b)) { bail!("invalid FASTA bases in {}",path.display()); }
            record.1.push_str(&seq);
        }
    }
    if records.is_empty() || records.iter().any(|r| r.1.is_empty()) { bail!("empty FASTA record in {}",path.display()); }
    let mut out=fs::File::create(output)?;
    for (i,(header,seq)) in records.iter().enumerate() {
        let name=format!("{sample}#1#contig_{}",i+1);
        mapping.push_str(&format!("{sample}\t{}\t{}\t{name}\n",path.display(),header.replace('\t'," ")));
        let mut tags=Vec::new();
        if header.split_whitespace().any(|v| matches!(v,"reads_decide=false"|"chosen_by=fewest_circular_molecules"|"unplaced=true"|"length_uncertain=true")) {
            tags.push("st:Z:unresolved_configuration".to_string());
        }
        if header.split_whitespace().any(|v| v.strip_prefix("reference_fill_bp=").is_some_and(|n| n.parse::<u64>().map_or(true, |x| x>0))) {
            tags.push("rf:Z:input_reference_fill".to_string());
        }
        if !tags.is_empty() { provenance.insert(name.clone(),tags); }

        writeln!(out,">{name}")?;
        for chunk in seq.as_bytes().chunks(80) { out.write_all(chunk)?;out.write_all(b"\n")?; }
    }
    Ok(())
}

// Prepared FASTA has one line per original sequence line. Compute the standard
// byte offsets, including newlines, for PGGB's mandatory input index.
fn fasta_index(path: &Path) -> Result<()> {
    let bytes=fs::read(path)?;
    let mut rows=String::new();
    let (mut name,mut length,mut offset,mut line_bases,mut line_bytes)=(String::new(),0usize,0usize,0usize,0usize);
    let mut cursor=0usize;
    for line in bytes.split_inclusive(|&b| b == b'\n') {
        if line.first()==Some(&b'>') {
            if !name.is_empty() { rows.push_str(&format!("{name}\t{length}\t{offset}\t{line_bases}\t{line_bytes}\n")); }
            name=String::from_utf8(line[1..].to_vec())?.trim().to_string();
            length=0;offset=cursor+line.len();line_bases=0;line_bytes=0;
        } else {
            let n=line.iter().filter(|&&b| b!=b'\n' && b!=b'\r').count();
            if line_bases==0 { line_bases=n;line_bytes=line.len(); }
            length+=n;
        }
        cursor+=line.len();
    }
    if !name.is_empty() { rows.push_str(&format!("{name}\t{length}\t{offset}\t{line_bases}\t{line_bytes}\n")); }
    fs::write(format!("{}.fai",path.display()),rows)?;Ok(())
}
fn gfas(dir: &Path, found: &mut Vec<PathBuf>) -> Result<()> {
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        if entry.file_type()?.is_dir() { gfas(&entry.path(), found)?; }
        else if entry.path().extension().is_some_and(|x| x == "gfa") { found.push(entry.path()); }
    }
    Ok(())
}
/// Exact generated PanSN names carry input evidence into normalized paths. If a
/// backend loses that correspondence, conservatively taint all its emitted paths.
fn propagate_provenance(path: &Path, provenance: &BTreeMap<String,Vec<String>>) -> Result<Value> {
    if provenance.is_empty() { return Ok(json!({"marked_paths":0,"fallback_all_paths":false})); }
    let text=fs::read_to_string(path)?;
    let names:BTreeSet<&str>=text.lines().filter_map(|l| { let f:Vec<_>=l.split('\t').collect(); (f.first()==Some(&"P") && f.len()>1).then(||f[1]) }).collect();
    let fallback=provenance.keys().any(|n| !names.contains(n.as_str()));
    if names.is_empty() { bail!("backend lost input provenance and emitted no paths to mark; refusing trusted normalization"); }
    let all:BTreeSet<String>=provenance.values().flatten().cloned().collect();
    let (mut out,mut count)=(String::new(),0usize);
    for line in text.lines() {
        let f:Vec<_>=line.split('\t').collect();
        out.push_str(line);
        if f.first()==Some(&"P") && f.len()>1 {
            let tags:Vec<&String>=if fallback {all.iter().collect()} else {provenance.get(f[1]).into_iter().flatten().collect()};
            if !tags.is_empty() { count+=1; }
            for tag in tags {
                let prefix=&tag[..3];
                // Preserve an existing same-key tag only if it is equally conservative.
                // Otherwise retain it separately and replace its interpretation safely.
                if !f.iter().any(|v| v.starts_with(prefix)) { out.push('\t');out.push_str(tag); }
                else if prefix=="st:" && !f.iter().any(|v| *v=="st:Z:unresolved_configuration") {
                    bail!("path {} has conflicting status provenance; refusing to overwrite or silently trust it",f[1]);
                }
                // Existing rf tags already identify reference-derived sequence: retain
                // their more specific source instead of replacing it with a generic label.
            }
        }
        out.push('\n');
    }
    fs::write(path,out)?;
    Ok(json!({"marked_paths":count,"fallback_all_paths":fallback,"reason":"input FASTA uncertainty/reference-fill annotations preserved; absent path correspondence taints all paths"}))
}
pub fn run(cfg: &Config) -> Result<Value> {
    if cfg.inputs.is_empty() || cfg.threads == 0 { bail!("provide at least one sample assembly and positive threads"); }
    if cfg.out.exists() { bail!("output already exists: {}", cfg.out.display()); }
    if !safe(&cfg.sample) { bail!("invalid graph sample name"); }
    let mut seen = BTreeSet::new();
    for (name,path) in &cfg.inputs {
        if !safe(name) || !seen.insert(name) { bail!("sample names must be unique PanSN-safe names: {name}"); }
        if !path.is_file() { bail!("input FASTA missing: {}", path.display()); }
    }
    let exe = cfg.executable.clone().unwrap_or_else(|| cfg.backend.name().into());
    let exe = if !exe.is_absolute() && exe.components().count() > 1 {
        fs::canonicalize(&exe).with_context(|| format!("resolve executable {}", exe.display()))?
    } else { exe };
    // An explicit environment's backend resolves sibling helpers locally, without
    // mutating the user's shell or global environment.
    let mut search: Vec<PathBuf> = std::env::split_paths(&std::env::var_os("PATH").unwrap_or_default()).collect();
    if let Some(parent) = exe.parent().filter(|p| !p.as_os_str().is_empty()) { search.insert(0,parent.to_path_buf()); }
    let path_env = std::env::join_paths(search)?;
    // Fail clearly before creating output when an optional backend is absent.
    let probe = Command::new(&exe).current_dir(fs::canonicalize(&cfg.inputs[0].1)?.parent().context("input has no parent directory")?).env("PATH",&path_env).arg(match cfg.backend { Backend::Pggb => "--version", Backend::Minigraph => "--version" })
        .output().with_context(|| format!("cannot execute {}: install the optional {} backend or specify its executable", exe.display(), cfg.backend.name()))?;
    if !probe.status.success() { bail!("{} version probe failed: {}", exe.display(), String::from_utf8_lossy(&probe.stderr)); }
    fs::create_dir(&cfg.out)?;
    let out = fs::canonicalize(&cfg.out)?;
    let external = out.join("external"); fs::create_dir(&external)?;
    let inputs = out.join("inputs"); fs::create_dir(&inputs)?;
    let mut mapping = String::from("sample\tinput\toriginal_header\tpath\n");
    let mut prepared = Vec::new();
    let mut provenance=BTreeMap::new();
    for (i,(sample,path)) in cfg.inputs.iter().enumerate() {
        let dest = inputs.join(format!("sample_{i}.fa")); prepare(path,sample,&dest,&mut mapping,&mut provenance)?; prepared.push(dest);
    }
    fs::write(out.join("input-names.tsv"),mapping)?;
    fs::write(out.join("backend-version.txt"),[probe.stdout,probe.stderr].concat())?;
    let raw = external.join("graph.gfa");
    let mut command = Command::new(&exe);
    command.current_dir(&out).env("PATH",&path_env);
    match cfg.backend {
        Backend::Minigraph => {
            command.arg("-cxggs").arg(format!("-t{}",cfg.threads)).args(&prepared)
                .stdout(fs::File::create(&raw)?);
        },
        Backend::Pggb => {
            let combined = inputs.join("all.fa"); let mut file = fs::File::create(&combined)?;
            for path in &prepared { std::io::copy(&mut fs::File::open(path)?, &mut file)?; }
            file.flush()?; fasta_index(&combined)?;
            command.arg("-i").arg(&combined).arg("-o").arg(&external).arg("-t").arg(cfg.threads.to_string())
                .arg("-n").arg(cfg.inputs.len().to_string()).stdout(fs::File::create(out.join("backend.stdout.log"))?);
        }
    }
    command.stdin(Stdio::null()).stderr(fs::File::create(out.join("backend.stderr.log"))?);
    let args: Vec<_> = command.get_args().map(|a|a.to_string_lossy().into_owned()).collect();
    fs::write(out.join("command.json"),serde_json::to_vec_pretty(&json!({"executable":exe,"args":args,"cwd":out,"inputs":cfg.inputs}))?)?;
    let status=command.status().context("launch external pan builder")?;
    fs::write(out.join("backend.exit"), status.code().map_or("signal".into(),|x|x.to_string()))?;
    if !status.success() { bail!("{} failed ({status}); see {}",cfg.backend.name(),out.join("backend.stderr.log").display()); }
    let graph = match cfg.backend {
        Backend::Minigraph => raw,
        Backend::Pggb => {
            let mut found=Vec::new(); gfas(&external,&mut found)?; found.sort();
            let final_graphs: Vec<_> = found.iter().filter(|p|p.file_name().unwrap().to_string_lossy().ends_with(".smooth.final.gfa")).cloned().collect();
            if final_graphs.len()==1 { final_graphs[0].clone() }
            else if final_graphs.is_empty() && found.len()==1 { found[0].clone() }
            else { bail!("PGGB produced {} GFA files and {} final graphs; refusing ambiguous output selection",found.len(),final_graphs.len()); }
        }
    };
    if fs::metadata(&graph)?.len()==0 { bail!("{} produced an empty graph",cfg.backend.name()); }
    let imported=crate::pan_external::run(&crate::pan_external::Config { input:graph.clone(),extra_paths:None,out:out.join("normalized"),sample:cfg.sample.clone() })?;
    let normalized=out.join("normalized/normalized.gfa");
    let propagated=propagate_provenance(&normalized,&provenance)?;
    fs::write(out.join("input-provenance.json"),serde_json::to_vec_pretty(&provenance)?)?;
    let report=json!({"provenance":propagated,"backend":cfg.backend.name(),"raw_graph":graph,"normalized_graph":out.join("normalized/normalized.gfa"),"input_names":out.join("input-names.tsv"),"inputs":cfg.inputs,"import":imported,"limitations":["Inputs must be assemblies, not raw WGS reads.","Minigraph does not recover all sample walks; rank-zero backbone is not evidence for other sample paths.","Backend parameters are default graph-shaping parameters; biological accuracy requires validation."]});
    fs::write(out.join("build.json"),serde_json::to_vec_pretty(&report)?)?; Ok(report)
}

#[cfg(all(test,unix))]
mod tests {
    use super::*;
    use std::os::unix::fs::PermissionsExt;
    fn fixture(label: &str, body:&str)->(PathBuf,Config) {
        let dir=std::env::temp_dir().join(format!("ovasm-backend-{label}-{}",std::process::id()));fs::create_dir_all(&dir).unwrap();
        let input=dir.join("sample.fa");fs::write(&input,">original name\nACGTACGTACGTACGTACGTACGTACGTACGT\n").unwrap();
        let exe=dir.join("backend");fs::write(&exe,format!("#!/bin/sh\nif [ \"$1\" = --version ] || [ \"$1\" = -v ]; then echo fake-test; exit 0; fi\n{body}\n")).unwrap();fs::set_permissions(&exe,fs::Permissions::from_mode(0o755)).unwrap();
        let cfg=Config{backend:Backend::Minigraph,inputs:vec![("sample".into(),input)],out:dir.join("result"),threads:2,executable:Some(exe),sample:"panel".into()};(dir,cfg)
    }
    #[test] fn prepared_fasta_index_has_correct_offsets() {
        let (d,c)=fixture("index","exit 0");
        let input=&c.inputs[0].1;
        fs::write(input,format!(">one\n{}\n{}\n>two\nTTAA\n","A".repeat(90),"C".repeat(15))).unwrap();
        let prepared=d.join("prepared.fa");prepare(input,"sample",&prepared,&mut String::new(),&mut BTreeMap::new()).unwrap();fasta_index(&prepared).unwrap();
        let text=fs::read_to_string(format!("{}.fai",prepared.display())).unwrap();
        let bytes=fs::read(&prepared).unwrap();
        for (i,line) in text.lines().enumerate() {
            let f:Vec<_>=line.split('\t').collect();let offset=f[2].parse::<usize>().unwrap();
            assert_eq!(bytes[offset],if i==0 {b'A'} else {b'T'});
            assert_eq!(f[1],if i==0 {"105"} else {"4"});
            assert_eq!(f[3],if i==0 {"80"} else {"4"});
        }
        fs::remove_dir_all(d).unwrap();
    }
    #[test] fn pggb_preserves_uncertain_fasta_headers_in_paths() {
        let script="while [ $# -gt 0 ]; do case \"$1\" in -o) out=$2; shift 2;; *) shift;; esac; done\nprintf 'H\\tVN:Z:1.0\\nS\\t1\\tACGTACGTACGTACGTACGTACGTACGTACGT\\nP\\tsample#1#contig_1\\t1+\\t*\\n' > \"$out/test.smooth.final.gfa\"";
        let (d,mut c)=fixture("pggb-provenance",script);
        c.backend=Backend::Pggb;
        fs::write(&c.inputs[0].1,">source reads_decide=false reference_fill_bp=12 length_uncertain=true\nACGTACGTACGTACGTACGTACGTACGTACGT\n").unwrap();
        run(&c).unwrap();
        let text=fs::read_to_string(c.out.join("normalized/normalized.gfa")).unwrap();
        assert!(text.lines().any(|l| l.starts_with("P\tsample#1#contig_1\t") && l.contains("st:Z:unresolved_configuration") && l.contains("rf:Z:input_reference_fill")));
        fs::remove_dir_all(d).unwrap();
    }
    #[test] fn lost_sample_names_conservatively_mark_backbone() {
        let (d,_)=fixture("fallback","exit 0");let g=d.join("graph.gfa");
        fs::write(&g,"H\tVN:Z:1.0\nS\t1\tACGT\trf:Z:original_source\nP\tbackbone#1#c\t1+\t*\n").unwrap();
        let p=BTreeMap::from([("sample#1#c".into(),vec!["st:Z:unresolved_configuration".into()])]);
        let report=propagate_provenance(&g,&p).unwrap();assert_eq!(report["fallback_all_paths"],true);
        let text=fs::read_to_string(g).unwrap();assert!(text.contains("rf:Z:original_source"));assert!(text.contains("st:Z:unresolved_configuration"));
        fs::remove_dir_all(d).unwrap();
    }
    #[test] fn backend_failure_is_not_success() { let (d,c)=fixture("failure","exit 7");let e=run(&c).unwrap_err().to_string();assert!(e.contains("failed"));assert!(!c.out.join("build.json").exists());fs::remove_dir_all(d).unwrap(); }
    #[test] fn empty_graph_is_not_success() { let (d,c)=fixture("empty","exit 0");assert!(run(&c).unwrap_err().to_string().contains("empty graph"));fs::remove_dir_all(d).unwrap(); }
    #[test] fn missing_backend_and_existing_output_fail_cleanly() { let (d,mut c)=fixture("missing","exit 0");c.executable=Some(d.join("missing"));assert!(run(&c).is_err());assert!(!c.out.exists());fs::create_dir(&c.out).unwrap();assert!(run(&c).unwrap_err().to_string().contains("already exists"));fs::remove_dir_all(d).unwrap(); }
    #[test] fn minigraph_dispatch_and_normalization() { let (d,c)=fixture("success","pwd > backend.cwd; printf 'H\\tVN:Z:1.0\\nS\\t1\\tACGTACGTACGTACGTACGTACGTACGTACGT\\nP\\tsample#1#contig_1\\t1+\\t*\\n'");let report=run(&c).unwrap();assert_eq!(report["backend"],"minigraph");assert!(c.out.join("normalized/normalized.gfa").exists());let command=fs::read_to_string(c.out.join("command.json")).unwrap();assert!(command.contains("-cxggs")&&command.contains("-t2"));assert_eq!(fs::read_to_string(c.out.join("backend.cwd")).unwrap().trim(),fs::canonicalize(&c.out).unwrap().to_str().unwrap());fs::remove_dir_all(d).unwrap(); }
}
