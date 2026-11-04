//! Strict GFA1/P/W and rGFA import. No sample paths are inferred from alignments.
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use anyhow::{bail, ensure, Context, Result};
use serde_json::{json, Value};

pub struct Config {
    pub input: PathBuf,
    pub extra_paths: Option<PathBuf>,
    pub out: PathBuf,
    pub sample: String,
}
#[derive(Clone)]
struct Segment { seq: String, tags: Vec<String> }
#[derive(Clone)]
struct Walk { name: String, steps: Vec<(String, bool)>, tags: Vec<String>, expected: Option<usize>, backbone: bool }
fn text(path: &Path) -> Result<String> {
    let f = fs::File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut r: Box<dyn Read> = if path.extension().is_some_and(|x| x == "gz") {
        Box::new(flate2::read::MultiGzDecoder::new(f))
    } else { Box::new(f) };
    let mut s=String::new(); r.read_to_string(&mut s)?; Ok(s)
}
fn safe(s:&str)->String {
    if crate::unify::valid_pansn_field(s) { return s.to_string(); }
    let mut out=String::from("encoded_");
    for b in s.bytes() { out.push_str(&format!("{b:02X}")); }
    out
}
fn name(s:&str, sample:&str)->String {
    let f:Vec<_>=s.split('#').collect();
    if f.len()==3 && f.iter().all(|x| !x.is_empty()) {
        format!("{}#{}#{}",safe(f[0]),safe(f[1]),safe(f[2]))
    } else { format!("{}#1#{}",safe(sample),safe(s)) }
}
fn orientation(s:&str)->Result<bool> { match s { "+"=>Ok(true), "-"=>Ok(false), _=>bail!("invalid orientation {s}") } }
fn psteps(s:&str)->Result<Vec<(String,bool)>> {
    ensure!(s!="*" && !s.contains(';'),"empty or jump-containing paths are unsupported");
    s.split(',').map(|x| { ensure!(x.len()>1,"empty path step"); let(n,o)=x.split_at(x.len()-1); Ok((n.into(),orientation(o)?)) }).collect()
}
fn wsteps(s:&str)->Result<Vec<(String,bool)>> {
    let mut out=Vec::new(); let b=s.as_bytes(); let mut i=0;
    while i<b.len() {
        ensure!(b[i]==b'>' || b[i]==b'<',"invalid W walk");
        let forward=b[i]==b'>'; i+=1; let start=i;
        while i<b.len() && b[i]!=b'>' && b[i]!=b'<' { i+=1; }
        ensure!(i>start,"empty W step"); out.push((s[start..i].into(),forward));
    }
    ensure!(!out.is_empty(),"empty W walk"); Ok(out)
}
fn parse_path(f:&[&str])->Result<Walk> {
    if f[0]=="P" {
        ensure!(f.len()>=4,"short P record"); let steps=psteps(f[2])?;
        if f[3]!="*" {
            let overlaps:Vec<_>=f[3].split(',').collect();
            ensure!(overlaps.len()==steps.len().saturating_sub(1) && overlaps.iter().all(|&x| x=="0M"),"P overlaps must be exactly 0M");
        }
        Ok(Walk{name:f[1].into(),steps,tags:f[4..].iter().map(|x|x.to_string()).collect(),expected:None,backbone:false})
    } else {
        ensure!(f.len()>=7,"short W record");
        ensure!(!f[1].is_empty() && !f[3].is_empty(), "missing W sample/sequence identity");
        let _:u64=f[2].parse().context("invalid W haplotype index")?;
        let expected=if f[4]!="*" && f[5]!="*" { Some(f[5].parse::<usize>()?.checked_sub(f[4].parse()?).context("reversed W coordinates")?) } else {None};
        // A W record is an interval, not necessarily a complete chromosome. Retain
        // interval identity explicitly; do not concatenate separate W records.
        let mut tags:Vec<String>=f[7..].iter().map(|x|x.to_string()).collect();
        tags.push(format!("wi:Z:{}:{}-{}",f[3],f[4],f[5]));
        if !tags.iter().any(|x|x.starts_with("st:")) { tags.push("st:Z:external_walk_interval".into()); }
        // A W interval is never proof of a complete sample path, irrespective
        // of any pre-existing status supplied by the external graph.
        if !tags.iter().any(|x|x.starts_with("up:")) { tags.push("up:Z:external_walk_interval".into()); }
        Ok(Walk{name:format!("{}#{}#{}:{}-{}",f[1],f[2],f[3],f[4],f[5]),steps:wsteps(f[6])?,tags,expected,backbone:false})
    }
}
fn tag<'a>(tags:&'a[String],prefix:&str)->Option<&'a str> { tags.iter().find_map(|x|x.strip_prefix(prefix)) }
fn normalize(input:&str, extra:Option<&str>, sample:&str)->Result<(String,Value,String)> {
    ensure!(crate::unify::valid_pansn_field(sample),"invalid sample label");
    let mut segs=BTreeMap::<String,Segment>::new();
    let mut links=Vec::<(String,bool,String,bool,Vec<String>)>::new();
    let mut walks=Vec::new();
    for line in input.lines().filter(|l| !l.is_empty() && !l.starts_with('#')) {
        let f:Vec<_>=line.split('\t').collect();
        match f[0] {
            "H"=>{},
            "S"=>{
                ensure!(f.len()>=3 && !f[1].is_empty(),"short S record");
                ensure!(f[2]!="*" && !f[2].is_empty() && f[2].bytes().all(|b|b"ACGTRYSWKMBDHVNacgtryswkmbdhvn".contains(&b)),"segment {} needs explicit IUPAC DNA sequence",f[1]);
                let tags:Vec<String>=f[3..].iter().map(|x|x.to_string()).collect();
                if let Some(n)=tag(&tags,"LN:i:") { ensure!(n.parse::<usize>()?==f[2].len(),"LN mismatch for {}",f[1]); }
                ensure!(segs.insert(f[1].into(),Segment{seq:f[2].into(),tags}).is_none(),"duplicate segment {}",f[1]);
            },
            "L"=>{
                ensure!(f.len()>=6 && f[5]=="0M","only explicit 0M links are supported; normalize overlaps first");
                links.push((f[1].into(),orientation(f[2])?,f[3].into(),orientation(f[4])?,f[6..].iter().map(|x|x.to_string()).collect()));
            },
            "P"|"W"=>walks.push(parse_path(&f)?),
            _=>bail!("unsupported GFA record {}; refusing to discard graph information",f[0]),
        }
    }
    if let Some(extra)=extra { for line in extra.lines().filter(|l|!l.is_empty() && !l.starts_with('#')) {
        let f:Vec<_>=line.split('\t').collect();
        ensure!(matches!(f[0],"P"|"W"),"additional paths file accepts only P/W records, never GAF alignments");
        walks.push(parse_path(&f)?);
    }}
    ensure!(!segs.is_empty(),"empty graph");
    let mut edges=BTreeSet::new();
    for (a,af,b,bf,_) in &links {
        ensure!(segs.contains_key(a) && segs.contains_key(b),"link references missing segment");
        edges.insert(((a.clone(),*af),(b.clone(),*bf)));
        edges.insert(((b.clone(),!*bf),(a.clone(),!*af)));
    }
    let supplied=walks.len();
    if walks.is_empty() {
        let mut groups=BTreeMap::<String,Vec<(usize,String)>>::new();
        for (id,s) in &segs {
            if tag(&s.tags,"SR:i:")==Some("0") {
                let sn=tag(&s.tags,"SN:Z:").context("rank-0 segment lacks SN")?;
                let so=tag(&s.tags,"SO:i:").context("rank-0 segment lacks SO")?.parse()?;
                groups.entry(sn.into()).or_default().push((so,id.clone()));
            }
        }
        ensure!(!groups.is_empty(),"no P/W paths or rGFA rank-0 backbone; provide explicit P/W paths");
        for (sn,mut ordered) in groups {
            ordered.sort(); let mut end=0usize; let mut steps=Vec::new();
            for (offset,id) in ordered {
                ensure!(offset==end,"rank-0 backbone {sn} has a gap/overlap or nonzero origin");
                end=offset.checked_add(segs[&id].seq.len()).context("backbone length overflow")?;
                steps.push((id,true));
            }
            walks.push(Walk{name:format!("{sample}#1#{sn}"),steps,tags:vec!["st:Z:reference_backbone_only".into()],expected:Some(end),backbone:true});
        }
    }
    let mut used=BTreeSet::new(); let mut names=BTreeSet::new();
    let mut mapping=String::from("original\tnormalized\tspelled_bases\trole\n");
    let mut out=format!("H\tVN:Z:1.0\tov:Z:organelleverse.graph.v1\tsn:Z:{sample}\n");
    for (id,s) in &segs { out.push_str(&format!("S\t{id}\t{}{}\n",s.seq,if s.tags.is_empty(){String::new()}else{format!("\t{}",s.tags.join("\t"))})); }
    let sign=|b|if b{"+"}else{"-"};
    for (a,af,b,bf,tags) in &links { out.push_str(&format!("L\t{a}\t{}\t{b}\t{}\t0M{}\n",sign(*af),sign(*bf),if tags.is_empty(){String::new()}else{format!("\t{}",tags.join("\t"))})); }
    for w in &walks {
        let normalized=name(&w.name,sample);
        ensure!(!segs.contains_key(&normalized), "path name {normalized} collides with a segment");
        ensure!(names.insert(normalized.clone()),"duplicate normalized path {normalized}; supply distinct interval/sample identities");
        let mut length=0usize;
        for (id,_) in &w.steps { let s=segs.get(id).with_context(||format!("path {} missing node {id}",w.name))?; length=length.checked_add(s.seq.len()).context("path length overflow")?; used.insert(id.clone()); }
        for pair in w.steps.windows(2) { ensure!(edges.contains(&(pair[0].clone(),pair[1].clone())),"path {} crosses absent link",w.name); }
        if w.tags.iter().any(|t|t=="TP:Z:circular") { ensure!(edges.contains(&(w.steps.last().unwrap().clone(),w.steps[0].clone())),"circular path {} lacks closing link",w.name); }
        if let Some(expected)=w.expected {ensure!(expected==length,"W/backbone length mismatch for {}",w.name);}
        let steps=w.steps.iter().map(|(id,f)|format!("{id}{}",sign(*f))).collect::<Vec<_>>().join(",");
        let tags=if w.tags.is_empty(){String::new()}else{format!("\t{}",w.tags.join("\t"))};
        out.push_str(&format!("P\t{normalized}\t{steps}\t*{tags}\n"));
        mapping.push_str(&format!("{}\t{normalized}\t{length}\t{}\n",w.name,if w.backbone{"reference_backbone"}else{"supplied_path"}));
    }
    let total:usize=segs.values().map(|s|s.seq.len()).sum();
    let covered:usize=used.iter().map(|id|segs[id].seq.len()).sum();
    let report=json!({"segments":segs.len(),"links":links.len(),"paths":walks.len(),"supplied_paths":supplied,
        "reference_backbone_only":supplied==0,"sample_path_coverage":{"covered_node_bases":covered,"stored_node_bases":total,"fraction":covered as f64/total as f64,"all_sample_paths_known":false},
        "sequence_and_path_lengths_validated":true,"physical_topology_certified":false,
        "limitations":["Import preserves graph paths, not evidence that reads resolve them.","rGFA rank-0 backbones do not represent paths of all samples; novel segments remain available for bait.","W records remain separate intervals; GAF partial alignments are not accepted as complete paths."]});
    Ok((out,report,mapping))
}
pub fn run(cfg:&Config)->Result<Value> {
    ensure!(!cfg.out.exists(),"output already exists: {}",cfg.out.display());
    let input=text(&cfg.input)?; let extra=cfg.extra_paths.as_deref().map(text).transpose()?;
    let (normalized,mut report,mapping)=normalize(&input,extra.as_deref(),&cfg.sample)?;
    if let Some(parent)=cfg.out.parent() {fs::create_dir_all(parent)?;}
    fs::create_dir(&cfg.out)?;
    let output=cfg.out.join("normalized.gfa");
    fs::write(&output,normalized)?; fs::write(cfg.out.join("path-names.tsv"),mapping)?;
    report["input"]=json!(cfg.input); report["extra_paths"]=json!(cfg.extra_paths); report["normalized_gfa"]=json!(output);
    fs::write(cfg.out.join("import.json"),serde_json::to_vec_pretty(&report)?)?; Ok(report)
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test] fn keeps_path_identity_tags_and_reverse_walks() {
        let g="S\t1\tACG\nS\t2\tTTA\nL\t1\t+\t2\t-\t0M\nP\ta#1#pt\t1+,2-\t0M\tst:Z:unresolved_configuration\n";
        let (out,r,_)=normalize(g,None,"panel").unwrap();
        assert!(out.contains("P\ta#1#pt\t1+,2-\t*\tst:Z:unresolved_configuration")); assert_eq!(r["supplied_paths"],1);
    }
    #[test] fn validates_w_interval_and_missing_edge() {
        let g="S\t1\tACG\nS\t2\tTTA\nL\t1\t+\t2\t-\t0M\nW\tsample\t1\tchr\t10\t16\t>1<2\n";
        let (out,_,_)=normalize(g,None,"panel").unwrap(); assert!(out.contains("external_walk_interval"));
        assert!(normalize(&g.replace("\t16\t", "\t17\t"),None,"panel").is_err());
        assert!(normalize(&g.replace("L\t1\t+\t2\t-\t0M\n",""),None,"panel").is_err());
    }
    #[test] fn rgfa_recovers_only_verified_backbone_and_keeps_novel_nodes() {
        let g="S\ta\tACG\tSN:Z:chr\tSO:i:0\tSR:i:0\nS\tb\tTTA\tSN:Z:chr\tSO:i:3\tSR:i:0\nS\tnovel\tCCC\tSR:i:1\nL\ta\t+\tb\t+\t0M\n";
        let (out,r,_)=normalize(g,None,"panel").unwrap(); assert!(out.contains("S\tnovel\tCCC")); assert_eq!(r["reference_backbone_only"],true);
        assert_eq!(r["sample_path_coverage"]["covered_node_bases"],6);
        assert!(normalize(&g.replace("SO:i:3","SO:i:4"),None,"panel").is_err());
    }
    #[test] fn rejects_missing_sequence_overlaps_and_pathless_graphs() {
        for g in ["S\t1\t*\n", "S\t1\tAAA\n", "S\t1\tAAA\nL\t1\t+\t1\t+\t1M\nP\tp\t1+\t*\n"] { assert!(normalize(g,None,"panel").is_err()); }
    }
    #[test] fn gzip_import_and_existing_output_refusal() {
        use std::io::Write;
        let dir=std::env::temp_dir().join(format!("ovasm-import-gzip-{}",std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let input=dir.join("input.gfa.gz");
        let mut z=flate2::write::GzEncoder::new(fs::File::create(&input).unwrap(),flate2::Compression::default());
        z.write_all(b"S\t1\tACGT\nP\tsample_one#1#pt\t1+\t*\trf:Z:source\n").unwrap(); z.finish().unwrap();
        let cfg=Config{input,extra_paths:None,out:dir.join("output"),sample:"container".into()};
        run(&cfg).unwrap();
        let original=fs::read(cfg.out.join("normalized.gfa")).unwrap();
        assert!(String::from_utf8_lossy(&original).contains("sample_one#1#pt"));
        assert!(run(&cfg).is_err());
        assert_eq!(original,fs::read(cfg.out.join("normalized.gfa")).unwrap());
        fs::remove_dir_all(dir).unwrap();
    }
    #[test] fn explicit_paths_are_checked_and_gaf_is_rejected() {
        let g="S\t1\tAAA\n";
        assert!(normalize(g,Some("P\tp\t1+\t*\n"),"panel").is_ok());
        assert!(normalize(g,Some("read\t100\t0\t90\t+\t>1\t100\t0\t90\t80\t90\t60"),"panel").is_err());
    }
}
