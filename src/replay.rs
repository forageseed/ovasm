//! Process-local, lossless replay of the original ordered recruitment batches.
//! Only sequences are cached. Final output replays the original input once for exact qualities.
//!
//! Caching the qualities too did not pay: deflated (8.5 GB for a 12.9 GB DToL HiFi library)
//! they cost the first pass about the 40 s they saved the final one, and stored as they are
//! (19.6 GB) the final pass took 335 s to read them from a busy disk, against 99 s to
//! decompress the input.

use std::fs::{File, OpenOptions};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use crate::recruit::{read_batches, Rec};
use anyhow::{ensure, Context, Result};
use rayon::prelude::*;

pub struct Entry {
    pub rec: Rec,
}

pub struct Replay {
    dir: PathBuf,
    sequence: Option<BufWriter<File>>,

    pub ready: bool,
    pub calls: Vec<crate::recruit::Call>,
    /// Per read, an upper bound on its coverage by any target since it was last classified
    /// (`recruit::still_other`).
    pub cover_bounds: Vec<u16>,
    pub reclassified_reads: Vec<u64>,
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct ReplayStats {
    pub source_classification_passes: usize,
    pub replay_classification_passes: usize,
    pub sequence_bytes: u64,
    pub final_output_input_passes: usize,
    pub reclassified_reads: Vec<u64>,
}

fn number(w: &mut impl Write, n: u64) -> Result<()> {
    w.write_all(&n.to_le_bytes())?;
    Ok(())
}
fn next_number(r: &mut impl Read) -> Result<u64> {
    let mut bytes = [0; 8];
    r.read_exact(&mut bytes)?;
    Ok(u64::from_le_bytes(bytes))
}

fn pack(seq: &[u8]) -> Option<Vec<u8>> {
    let mut out = vec![0; seq.len().div_ceil(4)];
    for (i, &b) in seq.iter().enumerate() {
        let value = match b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return None,
        };
        out[i / 4] |= value << (2 * (i % 4));
    }
    Some(out)
}

/// One read's `sequence.bin` record: id and sequence lengths, whether the sequence is packed,
/// the id, the sequence.
fn encode(r: &Rec) -> Vec<u8> {
    let packed = pack(&r.seq);
    let body = packed.as_deref().unwrap_or(&r.seq);
    let mut out = Vec::with_capacity(17 + r.id.len() + body.len());
    out.extend_from_slice(&(r.id.len() as u64).to_le_bytes());
    out.extend_from_slice(&(r.seq.len() as u64).to_le_bytes());
    out.push(u8::from(packed.is_some()));
    out.extend_from_slice(&r.id);
    out.extend_from_slice(body);
    out
}

/// One batch encoded for the cache (`Replay::encode`), written by `Replay::write`.
pub struct Encoded {
    records: Vec<Vec<u8>>,
}

/// The four bases of every packed byte (`pack` puts base i at bits 2(i % 4) of byte i / 4).
const UNPACKED: [[u8; 4]; 256] = {
    let mut table = [[0u8; 4]; 256];
    let mut byte = 0;
    while byte < 256 {
        let mut i = 0;
        while i < 4 {
            table[byte][i] = b"ACGT"[(byte >> (2 * i)) & 3];
            i += 1;
        }
        byte += 1;
    }
    table
};

/// Four bases per lookup: the cache is read and unpacked on one thread.
fn unpack(bytes: &[u8], len: usize) -> Vec<u8> {
    let mut seq = Vec::with_capacity(bytes.len() * 4);
    for &b in bytes {
        seq.extend_from_slice(&UNPACKED[b as usize]);
    }
    seq.truncate(len);
    seq
}

impl Replay {
    pub fn new(dir: &Path) -> Result<Self> {
        // No persistent-cache reuse: this directory belongs to one invocation.
        std::fs::create_dir(dir)
            .with_context(|| format!("recruit cache must be a new directory: {}", dir.display()))?;
        let create = |name| -> Result<BufWriter<File>> {
            Ok(BufWriter::with_capacity(
                1 << 20,
                OpenOptions::new()
                    .write(true)
                    .create_new(true)
                    .open(dir.join(name))?,
            ))
        };
        Ok(Self {
            dir: dir.into(),
            sequence: Some(create("sequence.bin")?),
            ready: false,
            calls: Vec::new(),
            cover_bounds: Vec::new(),
            reclassified_reads: Vec::new(),
        })
    }

    /// One batch for `write`, encoded on the worker threads: the scan encodes a batch while it
    /// classifies it (packing on the scan's own thread took 46-53 s of a ~200 s first pass).
    pub fn encode(batch: &[Entry]) -> Encoded {
        Encoded {
            records: batch.par_iter().map(|e| encode(&e.rec)).collect(),
        }
    }

    pub fn write(&mut self, batch: Encoded) -> Result<()> {
        ensure!(!self.ready, "recruit cache is already complete");
        let s = self.sequence.as_mut().expect("writer open");
        number(s, batch.records.len() as u64)?;
        for r in &batch.records {
            s.write_all(r)?;
        }
        Ok(())
    }

    pub fn store(&mut self, batch: &[Entry]) -> Result<()> {
        self.write(Self::encode(batch))
    }

    pub fn finish(&mut self) -> Result<()> {
        number(self.sequence.as_mut().expect("writer open"), 0)?;
        self.sequence.as_mut().expect("writer open").flush()?;
        self.sequence.take();
        self.ready = true;
        Ok(())
    }

    pub fn batches(
        &self,
    ) -> (
        crossbeam_channel::Receiver<Vec<Entry>>,
        std::thread::JoinHandle<Result<()>>,
    ) {
        let path = self.dir.join("sequence.bin");
        let (tx, rx) = crossbeam_channel::bounded(8);
        let thread = std::thread::spawn(move || -> Result<()> {
            let mut r = BufReader::with_capacity(1 << 20, File::open(path)?);
            loop {
                let tw = std::time::Instant::now();
                let count = next_number(&mut r)?;
                if count == 0 {
                    break;
                }
                let mut batch = Vec::with_capacity(count as usize);
                for _ in 0..count {
                    let id_len = next_number(&mut r)? as usize;
                    let len = next_number(&mut r)? as usize;
                    let mut flags = [0; 1];
                    r.read_exact(&mut flags)?;
                    ensure!(
                        flags.iter().all(|&f| f <= 1),
                        "invalid recruit cache encoding"
                    );
                    let mut id = vec![0; id_len];
                    r.read_exact(&mut id)?;
                    let mut seq = vec![0; if flags[0] == 1 { len.div_ceil(4) } else { len }];
                    r.read_exact(&mut seq)?;
                    if flags[0] == 1 {
                        seq = unpack(&seq, len);
                    }
                    batch.push(Entry {
                        rec: Rec {
                            id,
                            seq,
                            qual: None,
                        },
                    });
                }
                crate::recruit::PROF_READER_WORK_NS.fetch_add(tw.elapsed().as_nanos() as u64, std::sync::atomic::Ordering::Relaxed);
                let ts = std::time::Instant::now();
                if tx.send(batch).is_err() {
                    return Ok(());
                }
                crate::recruit::PROF_READER_SEND_NS.fetch_add(ts.elapsed().as_nanos() as u64, std::sync::atomic::Ordering::Relaxed);
            }
            Ok(())
        });
        (rx, thread)
    }

    /// The cached reads as plain records without qualities, as the depth sketch reads them.
    pub fn sequences(&self) -> crate::seedfree::Batches {
        let (cached, reader) = self.batches();
        let (tx, rx) = crossbeam_channel::bounded(2);
        let thread = std::thread::spawn(move || -> Result<()> {
            for batch in cached {
                if tx
                    .send(batch.into_iter().map(|entry| entry.rec).collect())
                    .is_err()
                {
                    break; // the consumer stopped; dropping `cached` stops the reader
                }
            }
            reader.join().expect("cache reader panicked")
        });
        (rx, thread)
    }

    pub fn stats(&self, passes: usize) -> Result<ReplayStats> {
        Ok(ReplayStats {
            source_classification_passes: 1,
            replay_classification_passes: passes - 1,
            sequence_bytes: std::fs::metadata(self.dir.join("sequence.bin"))?.len(),
            final_output_input_passes: usize::from(passes > 1),
            reclassified_reads: self.reclassified_reads.clone(),
        })
    }

    /// Called only after successful final output/report writing. Failed caches remain inspectable.
    pub fn remove(self) -> Result<()> {
        std::fs::remove_file(self.dir.join("sequence.bin"))?;
        std::fs::remove_dir(self.dir)?;
        Ok(())
    }
}

pub fn source_batches(
    inputs: &[PathBuf],
) -> (
    crossbeam_channel::Receiver<Vec<Entry>>,
    std::thread::JoinHandle<Result<()>>,
) {
    let (original, reader) = read_batches(inputs);
    let (tx, rx) = crossbeam_channel::bounded(1);
    let thread = std::thread::spawn(move || -> Result<()> {
        for batch in original {
            let batch = batch.into_iter().map(|rec| Entry { rec }).collect();
            if tx.send(batch).is_err() {
                return Ok(());
            }
        }
        reader.join().expect("input reader panicked")?;
        Ok(())
    });
    (rx, thread)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn cached_sequences_and_ids_are_exact_including_non_acgt_and_fasta() {
        let dir = std::env::temp_dir().join(format!("ovasm-replay-{}", std::process::id()));
        let mut cache = Replay::new(&dir).unwrap();
        let sequences: &[&[u8]] = &[b"ACGTA", b"acgtnNRY-", b"TT", b""];
        let entries: Vec<_> = sequences
            .iter()
            .enumerate()
            .map(|(i, seq)| Entry {
                rec: Rec {
                    id: format!("same name {i}").into_bytes(),
                    seq: seq.to_vec(),
                    qual: (i != 2).then(|| vec![b'!' + i as u8; seq.len()]),
                },
            })
            .collect();
        cache.store(&entries[..2]).unwrap();
        cache.store(&entries[2..]).unwrap();
        cache.finish().unwrap();
        let (rx, thread) = cache.sequences();
        let plain: Vec<Rec> = rx.iter().flatten().collect();
        thread.join().unwrap().unwrap();
        assert_eq!(plain.len(), entries.len());
        for (r, e) in plain.iter().zip(&entries) {
            assert_eq!((&r.id, &r.seq, &r.qual), (&e.rec.id, &e.rec.seq, &None));
        }
        let (rx, thread) = cache.batches();
        let mut index = 0;
        for batch in rx {
            for entry in batch {
                assert_eq!(entry.rec.id, entries[index].rec.id);
                assert_eq!(entry.rec.seq, entries[index].rec.seq);
                assert!(entry.rec.qual.is_none());
                index += 1;
            }
        }
        thread.join().unwrap().unwrap();
        assert_eq!(index, entries.len());
        cache.remove().unwrap();
        assert!(!dir.exists());
    }
    #[test]
    fn packing_round_trips_every_length() {
        for len in 0..40 {
            let seq: Vec<u8> = (0..len).map(|i| b"ACGT"[(i * 7 + len / 3) % 4]).collect();
            assert_eq!(unpack(&pack(&seq).unwrap(), len), seq);
        }
    }

    #[test]
    fn truncated_cache_fails_without_reopening_original_input() {
        let dir =
            std::env::temp_dir().join(format!("ovasm-replay-truncated-{}", std::process::id()));
        let mut cache = Replay::new(&dir).unwrap();
        cache.finish().unwrap();
        File::options()
            .write(true)
            .open(dir.join("sequence.bin"))
            .unwrap()
            .set_len(3)
            .unwrap();
        let (rx, thread) = cache.batches();
        assert!(rx.recv().is_err());
        assert!(thread.join().unwrap().is_err());
        cache.remove().unwrap();
    }
}
