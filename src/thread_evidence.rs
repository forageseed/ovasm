//! Single-record linkage for candidate scoring. One accepted complete record contributes
//! one maximal continuous walk; split chains, ambiguous assignments and non-ACGT spans
//! remain local evidence only. This is linkage evidence, not proof of molecular origin.
use std::collections::{BTreeMap, BTreeSet};
use std::hash::Hasher;

use anyhow::{bail, Result};
use rustc_hash::{FxHashMap, FxHasher};
use serde::Serialize;

use crate::evidence::{flip, Oriented, Step, Traversal};

pub const SCHEMA: &str = "ovasm.single-record-threads.v1";

#[derive(Serialize)]
pub struct ThreadLink {
    pub from: String,
    pub to: String,
    pub records: u64,
}

#[derive(Serialize)]
pub struct ThreadEvidence {
    pub schema: &'static str,
    pub inputs: Vec<String>,
    pub graph_fingerprint: String,
    pub record_identity: &'static str,
    pub contributing_records: u64,
    pub accepted_records: u64,
    pub duplicate_records: u64,
    pub rejected_ambiguous_records: u64,
    pub rejected_non_acgt_records: u64,
    pub rejected_disjoint_records: u64,
    /// Independently deduplicated local adjacency support; never implies repeat phasing.
    pub links: Vec<ThreadLink>,
    pub paths: Vec<Traversal>,
}

#[derive(Default)]
pub struct ThreadCounter {
    // Only contributing records enter the ledger. Memory is proportional to their number;
    // the complete sequence is fingerprinted, not retained. Short-read mates may share IDs.
    seen: FxHashMap<(Vec<u8>, Option<u64>), (usize, u64)>,
    counts: BTreeMap<Vec<Oriented>, u64>,
    links: BTreeMap<(Oriented, Oriented), u64>,
    accepted: u64,
    duplicates: u64,
    ambiguous: u64,
    non_acgt: u64,
    disjoint: u64,
}

fn canonical(w: &[Oriented]) -> Vec<Oriented> {
    let rc: Vec<_> = w.iter().rev().map(|&o| flip(o)).collect();
    if w <= rc.as_slice() { w.to_vec() } else { rc }
}

pub(crate) fn acgt_span(seq: &[u8], steps: &[Step]) -> bool {
    let Some(lo) = steps.iter().map(|s| s.2.min(s.3)).min() else { return false; };
    let hi = steps.iter().map(|s| s.2.max(s.3)).max().unwrap();
    let lo = (lo.max(0) as usize).min(seq.len());
    let hi = (hi.max(0) as usize).min(seq.len());
    hi > lo && seq[lo..hi].iter().all(|b| matches!(b.to_ascii_uppercase(), b'A'|b'C'|b'G'|b'T'))
}

impl ThreadCounter {
    pub fn observe(&mut self, id: &[u8], seq: &[u8], steps: &[Step], ambiguous: bool,
                   walks: &[Vec<Oriented>], extra_adjacencies: &[(Oriented, Oriented)],
                   short_reads: bool) -> Result<()> {
        if walks.is_empty() && steps.is_empty() && extra_adjacencies.is_empty() { return Ok(()); }
        let id = id.split(|b| b.is_ascii_whitespace()).next().unwrap_or_default();
        if id.is_empty() { bail!("single-record threads require nonempty original read IDs"); }
        let mut hasher = FxHasher::default();
        for b in seq { hasher.write_u8(b.to_ascii_uppercase()); }
        let signature = (seq.len(), hasher.finish());
        let key = (id.to_vec(), short_reads.then_some(signature.1));
        if let Some(previous) = self.seen.get(&key) {
            if *previous != signature {
                bail!("the same contributing read ID names different complete records in thread evidence");
            }
            self.duplicates += 1;
            return Ok(());
        }
        self.seen.insert(key, signature);
        let mut local = BTreeSet::new();
        for (a, b) in steps.iter().map(|s| (s.0, s.1)).chain(extra_adjacencies.iter().copied()) {
            local.insert((a,b).min((flip(b),flip(a))));
        }
        for edge in local { *self.links.entry(edge).or_default() += 1; }
        if walks.is_empty() { return Ok(()); }
        if ambiguous { self.ambiguous += 1; return Ok(()); }
        if !acgt_span(seq, steps) { self.non_acgt += 1; return Ok(()); }
        if walks.len() != 1 { self.disjoint += 1; return Ok(()); }
        if walks[0].len() < 3 { return Ok(()); }
        *self.counts.entry(canonical(&walks[0])).or_default() += 1;
        self.accepted += 1;
        Ok(())
    }

    pub fn report(self, graph: &crate::evidence::Graph, short_reads: bool, inputs: &[std::path::PathBuf]) -> ThreadEvidence {
        ThreadEvidence {
            schema: SCHEMA,
            inputs: inputs.iter().map(|p|p.display().to_string()).collect(),
            graph_fingerprint: crate::unzip::graph_fingerprint(graph),
            record_identity: if short_reads {
                "contributing original ID plus complete sequence fingerprint; same-name mates may be distinct records"
            } else {
                "contributing original ID; duplicate complete records count once; conflicting contributing IDs fail"
            },
            contributing_records: self.seen.len() as u64,
            accepted_records: self.accepted,
            duplicate_records: self.duplicates,
            rejected_ambiguous_records: self.ambiguous,
            rejected_non_acgt_records: self.non_acgt,
            rejected_disjoint_records: self.disjoint,
            links: self.links.into_iter().map(|((a,b),records)| ThreadLink {
                from: graph.label(a), to: graph.label(b), records,
            }).collect(),
            paths: self.counts.into_iter().map(|(w, reads)| Traversal {
                path: w.iter().map(|&o| graph.label(o)).collect::<Vec<_>>().join(" "), reads,
            }).collect(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn walk() -> Vec<Oriented> { vec![(0,true),(1,true),(2,true)] }
    fn steps() -> Vec<Step> { vec![((0,true),(1,true),1,4),((1,true),(2,true),4,7)] }
    fn add(c: &mut ThreadCounter, id: &[u8], seq: &[u8]) -> Result<()> {
        c.observe(id, seq, &steps(), false, &[walk()], &[], false)
    }
    #[test]
    fn duplicate_ids_cannot_inflate_support() {
        let mut c=ThreadCounter::default();
        add(&mut c,b"read1 old comment",b"ACGTACGT").unwrap();
        add(&mut c,b"read1 other comment",b"acgtacgt").unwrap();
        assert_eq!((c.accepted,c.duplicates),(1,1));
        assert_eq!(c.counts[&walk()],1);
    }
    #[test]
    fn conflicting_complete_records_fail() {
        let mut c=ThreadCounter::default();
        add(&mut c,b"read1",b"ACGTACGT").unwrap();
        assert!(add(&mut c,b"read1",b"ACGTACGA").is_err());
    }
    #[test]
    fn different_original_records_can_show_the_same_sequence() {
        let mut c=ThreadCounter::default();
        add(&mut c,b"read1",b"ACGTACGT").unwrap(); add(&mut c,b"read2",b"ACGTACGT").unwrap();
        assert_eq!(c.counts[&walk()],2);
    }
    #[test]
    fn unresolved_competition_is_not_a_joint_vote() {
        let mut c=ThreadCounter::default();
        c.observe(b"r",b"ACGTACGT",&steps(),true,&[walk()],&[],false).unwrap();
        assert_eq!((c.accepted,c.ambiguous),(0,1));
    }
    #[test]
    fn n_inside_crossing_span_is_excluded_but_not_outside() {
        let mut c=ThreadCounter::default();
        add(&mut c,b"bad",b"ACGNACGT").unwrap();
        add(&mut c,b"good",b"NCGTACGN").unwrap();
        assert_eq!((c.accepted,c.non_acgt),(1,1));
    }
    #[test]
    fn separated_chains_are_not_fused_into_a_record_path() {
        let mut c=ThreadCounter::default();
        c.observe(b"r",b"ACGTACGT",&steps(),false,&[walk(),walk()],&[],false).unwrap();
        assert_eq!((c.accepted,c.disjoint),(0,1));
    }
    #[test]
    fn reverse_complement_votes_for_the_same_linkage() {
        let mut c=ThreadCounter::default(); add(&mut c,b"a",b"ACGTACGT").unwrap();
        let rc=walk().iter().rev().map(|&o|flip(o)).collect();
        c.observe(b"b",b"ACGTACGT",&steps(),false,&[rc],&[],false).unwrap();
        assert_eq!(c.counts.len(),1); assert_eq!(c.counts[&walk()],2);
    }
    #[test]
    fn short_read_mates_with_shared_ids_remain_distinct() {
        let mut c=ThreadCounter::default();
        for seq in [b"ACGTACGT",b"ACGTACGA",b"ACGTACGT"] {
            c.observe(b"pair",seq,&steps(),false,&[walk()],&[],true).unwrap();
        }
        assert_eq!((c.accepted,c.duplicates),(2,1));
    }
    #[test]
    fn no_linkage_does_not_allocate_an_identity_ledger() {
        let mut c=ThreadCounter::default();
        c.observe(b"r",b"ACGTACGT",&[],false,&[],&[],false).unwrap();
        assert!(c.seen.is_empty()); assert_eq!(c.accepted,0);
    }
    #[test]
    fn repeated_local_edges_and_duplicate_records_do_not_inflate_support() {
        let mut c=ThreadCounter::default();
        let mut s=steps(); s.extend(steps());
        for _ in 0..2 { c.observe(b"r",b"ACGTACGT",&s,false,&[],&[],false).unwrap(); }
        assert_eq!(c.links.len(),2); assert!(c.links.values().all(|&n|n==1));
        assert_eq!((c.accepted,c.duplicates),(0,1));
    }
    #[test]
    fn exact_short_read_adjacencies_never_become_joint_evidence() {
        let mut c=ThreadCounter::default(); let edge=((0,true),(1,true));
        c.observe(b"r",b"ACGTACGT",&[],false,&[],&[edge],true).unwrap();
        assert_eq!(c.links[&edge],1); assert_eq!(c.accepted,0); assert!(c.counts.is_empty());
    }
}
