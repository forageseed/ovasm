//! Process-local, lossless replay of the original ordered recruitment batches.
//! Only sequences are cached. Final output replays the original input once for exact qualities.

use std::fs::{File, OpenOptions};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use crate::recruit::{read_batches, Rec};
use anyhow::{ensure, Context, Result};

pub struct Entry {
    pub rec: Rec,
}

pub struct Replay {
    dir: PathBuf,
    sequence: Option<BufWriter<File>>,

    pub ready: bool,
    pub calls: Vec<crate::recruit::Call>,
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

fn unpack(bytes: &[u8], len: usize) -> Vec<u8> {
    (0..len)
        .map(|i| b"ACGT"[((bytes[i / 4] >> (2 * (i % 4))) & 3) as usize])
        .collect()
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
            reclassified_reads: Vec::new(),
        })
    }

    pub fn store(&mut self, batch: &[Entry]) -> Result<()> {
        ensure!(!self.ready, "recruit cache is already complete");
        let s = self.sequence.as_mut().expect("writer open");
        number(s, batch.len() as u64)?;
        for entry in batch {
            let r = &entry.rec;
            let packed = pack(&r.seq);
            number(s, r.id.len() as u64)?;
            number(s, r.seq.len() as u64)?;
            s.write_all(&[u8::from(packed.is_some())])?;
            s.write_all(&r.id)?;
            s.write_all(packed.as_deref().unwrap_or(&r.seq))?;
        }
        Ok(())
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
                if tx.send(batch).is_err() {
                    return Ok(());
                }
            }
            Ok(())
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
