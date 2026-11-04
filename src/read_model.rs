//! Bounded, topology-neutral fit of candidate covers to graph-bound record evidence.
//! The lexicographic objective is explicit, not a calibrated likelihood or abundance.
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};

use anyhow::{bail, Context, Result};
use serde::Serialize;
use serde_json::Value;

use crate::evidence::{flip, Graph, Oriented};

pub const WORK_BUDGET: u64 = 10_000_000;
pub const MODEL: &str = "single_record_linkage_v1";

type Edge = (Oriented, Oriented);
fn edge(a: Oriented, b: Oriented) -> Edge { (a,b).min((flip(b),flip(a))) }

#[derive(Clone, Debug, Serialize, Default)]
pub struct ReadFit {
    pub joint_records_explained: u64,
    pub thread_records_explained: u64,
    pub unexplained_thread_records: u64,
    /// Each adjacency counts at most once per record and once per candidate cover.
    /// A record can support several adjacencies; this is NOT independent sample size.
    pub adjacency_record_occurrences_explained: u64,
    pub unobserved_junctions: usize,
    pub weak_junctions: usize,
    #[serde(skip)]
    pub matched_paths: Vec<usize>,
}

impl ReadFit {
    /// Less means a better evidence fit. Topology, length and molecule count do not enter.
    pub fn order(&self, other: &Self) -> Ordering {
        other.joint_records_explained.cmp(&self.joint_records_explained)
            .then(other.thread_records_explained.cmp(&self.thread_records_explained))
            .then(other.adjacency_record_occurrences_explained.cmp(&self.adjacency_record_occurrences_explained))
            .then(self.unobserved_junctions.cmp(&other.unobserved_junctions))
            .then(self.weak_junctions.cmp(&other.weak_junctions))
    }
    pub fn complete(&self) -> bool {
        self.unexplained_thread_records == 0 && self.unobserved_junctions == 0
            && self.weak_junctions == 0 && self.adjacency_record_occurrences_explained > 0
    }
}

pub struct ReadModel {
    paths: Vec<(Vec<Oriented>,u64,bool)>,
    links: BTreeMap<Edge,u64>,
    min_reads: u64,
}

#[derive(Clone, Debug, Serialize)]
pub struct ReadPathCover {
    pub scope: &'static str,
    /// Zero is the selected representative; other indices are in the full ranked list.
    pub candidate_indices: Vec<usize>,
    pub supported_thread_records: u64,
    pub represented_thread_records: u64,
    pub unrepresentable_thread_records: u64,
    pub low_support_thread_records: u64,
    pub minimum_family_size_proven: bool,
    pub budget_exhausted: bool,
}

fn exact_cover(masks: &[u64], target: u64, initial: Vec<usize>, budget: &mut u64) -> (Vec<usize>,bool) {
    struct Search<'a> {masks:&'a [u64],best:Vec<usize>,memo:BTreeMap<u64,usize>,exhausted:bool}
    fn visit(s:&mut Search,remaining:u64,path:&mut Vec<usize>,budget:&mut u64) {
        if !spend(budget,1) {s.exhausted=true;return;}
        if remaining==0 {if path.len()<s.best.len() {s.best=path.clone();}return;}
        if path.len()+1>=s.best.len() {return;}
        if s.memo.get(&remaining).is_some_and(|&depth|depth<=path.len()) {return;}
        s.memo.insert(remaining,path.len());
        if !spend(budget,s.masks.len()) {s.exhausted=true;return;}
        let gain=s.masks.iter().map(|m|(m&remaining).count_ones()).max().unwrap_or(0);
        if gain==0 || path.len()+(remaining.count_ones().div_ceil(gain) as usize)>=s.best.len() {return;}
        let pivot=1u64<<remaining.trailing_zeros();
        let mut choices:Vec<_>=s.masks.iter().enumerate().filter(|(_,m)|*m&pivot!=0)
            .map(|(i,m)|(i,(m&remaining).count_ones())).collect();
        choices.sort_by_key(|&(i,g)|(std::cmp::Reverse(g),i));
        for (i,_) in choices {
            path.push(i);visit(s,remaining&!s.masks[i],path,budget);path.pop();
            if s.exhausted {return;}
        }
    }
    let mut search=Search{masks,best:initial,memo:BTreeMap::new(),exhausted:false};
    visit(&mut search,target,&mut Vec::new(),budget);
    (search.best,!search.exhausted)
}

fn spend(budget: &mut u64, cost: usize) -> bool {
    match budget.checked_sub(cost as u64) { Some(left)=>{*budget=left;true},None=>false }
}

fn contains(path: &[Oriented], circular: bool, query: &[Oriented], budget: &mut u64) -> Option<bool> {
    if path.is_empty() || query.is_empty() { return Some(false); }
    if !circular && query.len()>path.len() { return Some(false); }
    let rc: Vec<_> = query.iter().rev().map(|&x|flip(x)).collect();
    for q in [query,rc.as_slice()] {
        let starts=if circular {path.len()} else {path.len()-q.len()+1};
        for start in 0..starts {
            let mut matches=true;
            for (i,&o) in q.iter().enumerate() {
                if !spend(budget,1) { return None; }
                if path[(start+i)%path.len()]!=o { matches=false;break; }
            }
            if matches { return Some(true); }
        }
    }
    Some(false)
}

impl ReadModel {
    pub fn load(g: &Graph, report: &Value, repeats: &[bool], min_reads: u64) -> Result<Option<Self>> {
        let Some(v)=report.get("thread_evidence") else { return Ok(None); };
        if v["schema"].as_str()!=Some(crate::thread_evidence::SCHEMA) {
            bail!("unknown single-record evidence schema; regenerate check evidence");
        }
        if v["graph_fingerprint"].as_str()!=Some(crate::unzip::graph_fingerprint(g).as_str()) {
            bail!("single-record evidence does not match graph sequences, links and overlaps");
        }
        let index: BTreeMap<_,_>=g.names.iter().enumerate().map(|(i,n)|(n.as_str(),i)).collect();
        let parse=|s:&str| -> Result<Oriented> {
            let Some((name,sign))=s.len().checked_sub(1).map(|n|s.split_at(n)) else {bail!("empty thread label")};
            if sign!="+" && sign!="-" {bail!("invalid thread orientation")}
            Ok((*index.get(name).context("thread names unknown graph segment")?,sign=="+"))
        };
        let valid: BTreeSet<_>=g.links.iter().map(|&(a,b)|edge(a,b)).collect();
        let mut links=BTreeMap::new();
        for l in v["links"].as_array().context("thread evidence has no local links")? {
            let e=edge(parse(l["from"].as_str().context("missing link from")?)?,parse(l["to"].as_str().context("missing link to")?)?);
            let n=l["records"].as_u64().context("invalid thread link count")?;
            if !valid.contains(&e) || links.insert(e,n).is_some() {bail!("unknown or duplicate thread evidence link")}
        }
        let mut paths=Vec::new(); let mut seen=BTreeSet::new(); let mut total=0u64;
        for p in v["paths"].as_array().context("thread evidence has no paths")? {
            let w=p["path"].as_str().context("missing thread path")?.split_whitespace().map(parse).collect::<Result<Vec<_>>>()?;
            let n=p["reads"].as_u64().filter(|&n|n>0).context("invalid thread record count")?;
            if w.len()<3 || w.windows(2).any(|x|!valid.contains(&edge(x[0],x[1]))) {bail!("thread path does not walk the graph")}
            let rc:Vec<_>=w.iter().rev().map(|&x|flip(x)).collect();
            let key=w.clone().min(rc);
            if !seen.insert(key) {bail!("duplicate canonical record thread")}
            // Two complete repeat visits can concern the same repeat copy label. End
            // nodes are only partial anchors and must not manufacture joint phasing.
            let joint=w[1..w.len()-1].iter().filter(|o|repeats[o.0]).count()>=2;
            total=total.checked_add(n).context("thread counts overflow")?;
            paths.push((w,n,joint));
        }
        if Some(total)!=v["accepted_records"].as_u64() {bail!("single-record thread totals disagree")}
        Ok(Some(Self{paths,links,min_reads:min_reads.max(1)}))
    }

    pub fn for_anchors(&self, anchors: &[usize], repeats: &[bool]) -> Self {
        let belongs:BTreeSet<_>=anchors.iter().copied().collect();
        Self {paths:self.paths.iter().filter(|(w,_,_)|
            w.iter().any(|o|!repeats[o.0] && belongs.contains(&o.0))).cloned().collect(),
            links:self.links.clone(),min_reads:self.min_reads}
    }

    /// A minimal explanatory family is a representation of observed paths, not a count
    /// or identification of physical isoforms. Exact cardinality search is bounded.
    pub fn path_cover(&self, fits: &[&ReadFit], mut budget: u64) -> ReadPathCover {
        let eligible:BTreeSet<_>=self.paths.iter().enumerate().filter(|(_,(_,n,_))|*n>=self.min_reads).map(|(i,_)|i).collect();
        let target:BTreeSet<_>=fits.iter().flat_map(|f|f.matched_paths.iter().copied()).filter(|i|eligible.contains(i)).collect();
        let mut remaining=target.clone();let mut selected=Vec::new();let mut exhausted=false;
        while !remaining.is_empty() {
            let mut best=None;
            for (i,fit) in fits.iter().enumerate() {
                if !spend(&mut budget,fit.matched_paths.len()+1) {exhausted=true;break;}
                let gain=fit.matched_paths.iter().filter(|i|remaining.contains(i)).count();
                if gain>0 && best.map_or(true,|(_,old)|gain>old) {best=Some((i,gain));}
            }
            if exhausted {break;}
            let Some((i,_))=best else {break;};
            selected.push(i);for j in &fits[i].matched_paths {remaining.remove(j);}
        }
        let mut proven=target.is_empty();
        if !target.is_empty() && remaining.is_empty() && target.len()<=64 && !exhausted {
            let bit:BTreeMap<_,_>=target.iter().enumerate().map(|(i,&j)|(j,i)).collect();
            let masks:Vec<u64>=fits.iter().map(|f|f.matched_paths.iter().filter_map(|j|bit.get(j)).fold(0u64,|mask,i|mask|(1u64<<i))).collect();
            let all=masks.iter().fold(0,|a,m|a|m);
            (selected,proven)=exact_cover(&masks,all,selected,&mut budget);
            exhausted=!proven;
        }
        selected.sort_unstable();
        let represented:BTreeSet<_>=selected.iter().flat_map(|&i|fits[i].matched_paths.iter().copied()).collect();
        let count=|indices:&BTreeSet<usize>|indices.iter().map(|&i|self.paths[i].1).sum::<u64>();
        ReadPathCover {scope:"minimum explanatory family within enumerated candidates; not physical isoform count or abundance",
            candidate_indices:selected,supported_thread_records:count(&eligible),represented_thread_records:count(&represented),
            unrepresentable_thread_records:count(&eligible)-count(&target),
            low_support_thread_records:self.paths.iter().filter(|(_,n,_)|*n<self.min_reads).map(|(_,n,_)|*n).sum(),
            minimum_family_size_proven:proven,budget_exhausted:exhausted}
    }

    /// None means budget exhausted: callers must discard ALL partial new-model ranks.
    pub fn score(&self, walks: &[Vec<Oriented>], circular: &[bool], budget: &mut u64) -> Option<ReadFit> {
        let mut fit=ReadFit::default();
        let mut used=BTreeSet::new();
        for (path,&closed) in walks.iter().zip(circular) {
            if !spend(budget,path.len()) {return None;}
            for pair in path.windows(2) {used.insert(edge(pair[0],pair[1]));}
            if closed && !path.is_empty() {used.insert(edge(*path.last().unwrap(),path[0]));}
        }
        for e in used {
            let n=self.links.get(&e).copied().unwrap_or(0);
            if n==0 {fit.unobserved_junctions+=1;}
            else if n<self.min_reads {fit.weak_junctions+=1;}
            else {fit.adjacency_record_occurrences_explained+=n;}
        }
        for (index,(w,n,joint)) in self.paths.iter().enumerate() {
            let mut matched=false;
            for (p,&c) in walks.iter().zip(circular) {
                if contains(p,c,w,budget)? {matched=true;break;}
            }
            if !matched {fit.unexplained_thread_records+=n;}
            else if *n>=self.min_reads {
                fit.matched_paths.push(index);
                fit.thread_records_explained+=n;
                if *joint {fit.joint_records_explained+=n;}
            }
        }
        Some(fit)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn w(ns:&[usize])->Vec<Oriented>{ns.iter().map(|&n|(n,true)).collect()}
    #[test]
    fn wrapping_requires_closed_path_and_works_on_both_strands() {
        let p=w(&[0,1,2]);let q=w(&[2,0,1]);let mut b=100;
        assert_eq!(contains(&p,false,&q,&mut b),Some(false));
        assert_eq!(contains(&p,true,&q,&mut b),Some(true));
        let rc:Vec<_>=q.iter().rev().map(|&o|flip(o)).collect();
        assert_eq!(contains(&p,true,&rc,&mut b),Some(true));
    }
    #[test]
    fn equal_local_marginals_do_not_substitute_for_joint_linkage() {
        let model=ReadModel{paths:vec![(w(&[0,1,2,3,4]),10,true),(w(&[5,1,2,3,6]),10,true)],links:BTreeMap::new(),min_reads:2};
        let mut budget=10000;
        let yes=model.score(&[w(&[0,1,2,3,4]),w(&[5,1,2,3,6])],&[false,false],&mut budget).unwrap();
        let no=model.score(&[w(&[0,1,2,3,6]),w(&[5,1,2,3,4])],&[false,false],&mut budget).unwrap();
        assert_eq!((yes.joint_records_explained,no.joint_records_explained),(20,0));
        assert_eq!(yes.order(&no),Ordering::Less);
    }
    #[test]
    fn coexisting_configurations_are_not_forced_apart() {
        let model=ReadModel{paths:vec![(w(&[0,1,2]),10,true),(w(&[0,1,3]),10,true)],links:BTreeMap::new(),min_reads:2};
        let mut budget=1000;
        let a=model.score(&[w(&[0,1,2])],&[false],&mut budget).unwrap();
        let b=model.score(&[w(&[0,1,3])],&[false],&mut budget).unwrap();
        assert_eq!(a.order(&b),Ordering::Equal);assert_eq!(a.unexplained_thread_records,10);assert!(!a.complete());
    }
    #[test]
    fn closure_has_no_bonus_without_observed_connection() {
        let model=ReadModel{paths:vec![],links:[(edge((0,true),(1,true)),10),(edge((1,true),(2,true)),10)].into_iter().collect(),min_reads:2};
        let mut budget=1000;
        let open=model.score(&[w(&[0,1,2])],&[false],&mut budget).unwrap();
        let closed=model.score(&[w(&[0,1,2])],&[true],&mut budget).unwrap();
        assert_eq!(open.order(&closed),Ordering::Less);assert_eq!(closed.unobserved_junctions,1);
    }
    #[test]
    fn closing_evidence_supplies_an_actual_gain() {
        let links=[(edge((0,true),(1,true)),10),(edge((1,true),(2,true)),10),(edge((2,true),(0,true)),10)].into_iter().collect();
        let model=ReadModel{paths:vec![],links,min_reads:2};let mut b=1000;
        let open=model.score(&[w(&[0,1,2])],&[false],&mut b).unwrap();
        let closed=model.score(&[w(&[0,1,2])],&[true],&mut b).unwrap();
        assert_eq!(closed.order(&open),Ordering::Less);
    }
    #[test]
    fn repeat_loops_do_not_recount_the_same_edge_or_record() {
        let model=ReadModel{paths:vec![(w(&[0,1,0]),10,false)],links:[(edge((0,true),(1,true)),10),(edge((1,true),(0,true)),10)].into_iter().collect(),min_reads:2};
        let mut b=1000;
        let a=model.score(&[w(&[0,1,0])],&[false],&mut b).unwrap();
        let z=model.score(&[w(&[0,1,0,1,0])],&[false],&mut b).unwrap();
        assert_eq!(a.thread_records_explained,10);assert_eq!(a.order(&z),Ordering::Equal);
    }
    #[test]
    fn one_record_cannot_become_supported_joint_structure() {
        let model=ReadModel{paths:vec![(w(&[0,1,2]),1,true)],links:BTreeMap::new(),min_reads:2};let mut b=1000;
        let a=model.score(&[w(&[0,1,2])],&[false],&mut b).unwrap();
        assert_eq!((a.joint_records_explained,a.thread_records_explained),(0,0));
    }
    #[test]
    fn exhausted_budget_never_returns_a_partial_score() {
        let model=ReadModel{paths:vec![(w(&[0,1,2]),10,false)],links:BTreeMap::new(),min_reads:2};let mut b=1;
        assert!(model.score(&[w(&[0,1,2])],&[false],&mut b).is_none());
    }
    #[test]
    fn mismatched_graph_is_rejected_before_scoring() {
        let g=Graph{names:vec!["A".into()],seqs:vec![b"ACGT".to_vec()],links:vec![],overlaps:vec![]};
        let report=serde_json::json!({"thread_evidence":{"schema":crate::thread_evidence::SCHEMA,"graph_fingerprint":"wrong"}});
        assert!(ReadModel::load(&g,&report,&[false],2).err().unwrap().to_string().contains("does not match graph"));
    }
    #[test]
    fn another_components_records_are_not_unexplained_evidence() {
        let model=ReadModel{paths:vec![(w(&[0,1,2]),10,false),(w(&[3,4,5]),20,false)],links:BTreeMap::new(),min_reads:2};
        let scoped=model.for_anchors(&[0,1,2],&[false;6]);let mut b=1000;
        let fit=scoped.score(&[w(&[0,1,2])],&[false],&mut b).unwrap();
        assert_eq!(fit.thread_records_explained,10);assert_eq!(fit.unexplained_thread_records,0);
    }
    #[test]
    fn explanatory_family_uses_exact_search_when_greedy_is_suboptimal() {
        let model=ReadModel{paths:(0..6).map(|i|(w(&[i]),10,false)).collect(),links:BTreeMap::new(),min_reads:2};
        let f=|xs:&[usize]|ReadFit{matched_paths:xs.to_vec(),..ReadFit::default()};
        let fits=[f(&[0,1,2,3]),f(&[0,1,4]),f(&[2,3,5]),f(&[4]),f(&[5])];
        let cover=model.path_cover(&fits.iter().collect::<Vec<_>>(),10000);
        assert_eq!(cover.candidate_indices,vec![1,2]);assert!(cover.minimum_family_size_proven);
        assert_eq!(cover.represented_thread_records,60);
    }
    #[test]
    fn explanatory_family_preserves_unrepresentable_and_weak_records() {
        let model=ReadModel{paths:vec![(w(&[0]),10,false),(w(&[1]),5,false),(w(&[2]),1,false)],links:BTreeMap::new(),min_reads:2};
        let fit=ReadFit{matched_paths:vec![0],..ReadFit::default()};let c=model.path_cover(&[&fit],1000);
        assert_eq!((c.supported_thread_records,c.represented_thread_records,c.unrepresentable_thread_records,c.low_support_thread_records),(15,10,5,1));
    }
    #[test]
    fn explanatory_family_budget_does_not_manufacture_optimality() {
        let model=ReadModel{paths:vec![(w(&[0]),10,false)],links:BTreeMap::new(),min_reads:2};
        let fit=ReadFit{matched_paths:vec![0],..ReadFit::default()};let c=model.path_cover(&[&fit],0);
        assert!(c.budget_exhausted);assert!(!c.minimum_family_size_proven);assert_eq!(c.represented_thread_records,0);
    }
}
