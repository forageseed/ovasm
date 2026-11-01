//! Optional CUDA offload of the replay passes' check for changed k-mers (`OVASM_GPU=1`, a
//! prototype). Every replay pass asks of every cached read whether it holds a k-mer the last
//! extension changed (`BaitSet::any_key`); that is most of a pass's CPU time once few reads are
//! reclassified, and it is a GPU's kind of work: billions of independent lookups in a table of
//! megabytes. The kernel (`gpu_delta.cu`) uses the same index bit for bit, so the answer, and with
//! it every output, is the CPU's.
//!
//! The CUDA driver is opened at run time (`dlopen("libcuda.so.1")`) and the kernel is shipped as
//! PTX the driver compiles for the GPU at hand: ovasm builds and runs without CUDA, and falls back
//! to the CPU when there is no driver or GPU (Linux only).

use anyhow::{bail, ensure, Context, Result};

use crate::bait::BaitSet;
use crate::replay::Entry;

const PTX: &str = concat!(include_str!("gpu_delta.ptx"), "\0");
/// K-mer starts per GPU thread.
const CHUNK: u32 = 256;
const BLOCK_THREADS: u32 = 256;

/// Whether `OVASM_GPU=1` asks for the GPU.
pub fn requested() -> bool {
    std::env::var("OVASM_GPU").is_ok_and(|v| v == "1")
}

#[cfg(target_os = "linux")]
mod driver {
    use std::ffi::c_void;

    pub type CuResult = i32;
    pub type DevPtr = u64;

    extern "C" {
        fn dlopen(file: *const i8, mode: i32) -> *mut c_void;
        fn dlsym(handle: *mut c_void, name: *const i8) -> *mut c_void;
    }
    const RTLD_NOW: i32 = 2;

    macro_rules! api {
        ($($field:ident = $sym:literal : fn($($arg:ty),*);)*) => {
            pub struct Api { $(pub $field: unsafe extern "C" fn($($arg),*) -> CuResult,)* }
            impl Api {
                pub fn open() -> Option<Self> {
                    // SAFETY: dlopen/dlsym on NUL-terminated names; each symbol is the driver
                    // API function of that name and signature (cuda.h, API version 2 names).
                    unsafe {
                        let lib = dlopen(c"libcuda.so.1".as_ptr(), RTLD_NOW);
                        if lib.is_null() {
                            return None;
                        }
                        Some(Self { $($field: {
                            let f = dlsym(lib, concat!($sym, "\0").as_ptr() as *const i8);
                            if f.is_null() { return None; }
                            std::mem::transmute::<*mut c_void, unsafe extern "C" fn($($arg),*) -> CuResult>(f)
                        },)* })
                    }
                }
            }
        };
    }

    api! {
        init = "cuInit": fn(u32);
        device_get = "cuDeviceGet": fn(*mut i32, i32);
        primary_ctx_retain = "cuDevicePrimaryCtxRetain": fn(*mut *mut c_void, i32);
        ctx_set_current = "cuCtxSetCurrent": fn(*mut c_void);
        module_load_data = "cuModuleLoadData": fn(*mut *mut c_void, *const c_void);
        module_get_function = "cuModuleGetFunction": fn(*mut *mut c_void, *mut c_void, *const i8);
        mem_alloc = "cuMemAlloc_v2": fn(*mut DevPtr, usize);
        mem_free = "cuMemFree_v2": fn(DevPtr);
        htod = "cuMemcpyHtoD_v2": fn(DevPtr, *const c_void, usize);
        dtoh = "cuMemcpyDtoH_v2": fn(*mut c_void, DevPtr, usize);
        memset_d8 = "cuMemsetD8_v2": fn(DevPtr, u8, usize);
        launch = "cuLaunchKernel": fn(*mut c_void, u32, u32, u32, u32, u32, u32, u32, *mut c_void, *mut *mut c_void, *mut *mut c_void);
        synchronize = "cuCtxSynchronize": fn();
    }
}

#[cfg(target_os = "linux")]
use driver::{Api, DevPtr};

/// A device buffer that grows to the largest request and is reused.
#[cfg(target_os = "linux")]
#[derive(Default)]
struct Buf {
    ptr: DevPtr,
    cap: usize,
}

#[cfg(target_os = "linux")]
pub struct Gpu {
    api: Api,
    ctx: *mut std::ffi::c_void,
    func: *mut std::ffi::c_void,
    /// packed, starts, lens, chunk_read, chunk_first, bloom, keys, hit
    bufs: [Buf; 8],
    bloom_mask: u64,
    table_mask: u64,
    k: u32,
}

#[cfg(target_os = "linux")]
fn check(r: driver::CuResult, what: &str) -> Result<()> {
    ensure!(r == 0, "CUDA {what} failed (error {r})");
    Ok(())
}

#[cfg(target_os = "linux")]
impl Gpu {
    /// The first GPU (`OVASM_GPU_DEVICE` picks another) with the kernel loaded.
    pub fn open() -> Result<Self> {
        let api = Api::open().context("no CUDA driver (libcuda.so.1)")?;
        let device: i32 = std::env::var("OVASM_GPU_DEVICE").ok().and_then(|v| v.parse().ok()).unwrap_or(0);
        // SAFETY: driver API calls with valid out-pointers; the PTX is NUL-terminated.
        unsafe {
            check((api.init)(0), "init")?;
            let mut dev = 0;
            check((api.device_get)(&mut dev, device), "device")?;
            let mut ctx = std::ptr::null_mut();
            check((api.primary_ctx_retain)(&mut ctx, dev), "context")?;
            check((api.ctx_set_current)(ctx), "context")?;
            let mut module = std::ptr::null_mut();
            check((api.module_load_data)(&mut module, PTX.as_ptr().cast()), "kernel load")?;
            let mut func = std::ptr::null_mut();
            check((api.module_get_function)(&mut func, module, c"any_delta".as_ptr()), "kernel lookup")?;
            Ok(Self { api, ctx, func, bufs: Default::default(), bloom_mask: 0, table_mask: 0, k: 0 })
        }
    }

    fn upload(&mut self, buf: usize, bytes: &[u8]) -> Result<DevPtr> {
        let b = &mut self.bufs[buf];
        // SAFETY: device allocations owned by `self`; host slices are valid for their length.
        unsafe {
            if b.cap < bytes.len().max(1) {
                if b.cap > 0 {
                    // forgotten before the allocation, which can fail: a later pass would
                    // otherwise copy into the freed address, perhaps another buffer's by then
                    b.cap = 0;
                    check((self.api.mem_free)(b.ptr), "free")?;
                }
                let cap = bytes.len().max(1).next_power_of_two();
                check((self.api.mem_alloc)(&mut b.ptr, cap), "alloc")?;
                b.cap = cap;
            }
            if !bytes.is_empty() {
                check((self.api.htod)(b.ptr, bytes.as_ptr().cast(), bytes.len()), "copy to GPU")?;
            }
        }
        Ok(b.ptr)
    }

    /// The change set a pass checks reads against (`BaitSet::take_changes`).
    pub fn set_delta(&mut self, delta: &BaitSet) -> Result<()> {
        // SAFETY: the context belongs to this process.
        unsafe { check((self.api.ctx_set_current)(self.ctx), "context")? };
        let Some((bloom, bloom_mask, keys, table_mask)) = delta.lookup_parts() else {
            bail!("change set has no prepared index");
        };
        self.upload(5, as_bytes(bloom))?;
        self.upload(6, as_bytes(keys))?;
        (self.bloom_mask, self.table_mask, self.k) = (bloom_mask, table_mask, delta.k as u32);
        Ok(())
    }

    /// Per entry, whether it holds a k-mer of the change set, for the reads cached packed;
    /// `None` for the others (non-ACGT reads), which the CPU checks.
    pub fn changed(&mut self, batch: &[Entry]) -> Result<Vec<Option<bool>>> {
        let k = self.k;
        let (mut packed, mut starts, mut lens, mut chunk_read, mut chunk_first) =
            (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
        let mut slot = vec![None; batch.len()];
        for (i, e) in batch.iter().enumerate() {
            let Some(p) = &e.packed else { continue };
            let len = e.rec.seq.len() as u32;
            if len < k {
                slot[i] = Some(u32::MAX);
                continue;
            }
            let r = starts.len() as u32;
            slot[i] = Some(r);
            starts.push(packed.len() as u64);
            lens.push(len);
            packed.extend_from_slice(p);
            for first in (0..len - k + 1).step_by(CHUNK as usize) {
                chunk_read.push(r);
                chunk_first.push(first);
            }
        }
        let n_reads = starts.len();
        let mut hit = vec![0u8; n_reads];
        if n_reads > 0 {
            // SAFETY: the context belongs to this process; kernel arguments point at locals of
            // the types `any_delta` takes, live until the synchronizing copy back.
            unsafe {
                check((self.api.ctx_set_current)(self.ctx), "context")?;
                let mut d_packed = self.upload(0, &packed)?;
                let mut d_starts = self.upload(1, as_bytes(&starts))?;
                let mut d_lens = self.upload(2, as_bytes(&lens))?;
                let mut d_chunk_read = self.upload(3, as_bytes(&chunk_read))?;
                let mut d_chunk_first = self.upload(4, as_bytes(&chunk_first))?;
                let mut d_hit = self.upload(7, &hit)?;
                check((self.api.memset_d8)(d_hit, 0, n_reads), "clear")?;
                let (mut d_bloom, mut d_keys) = (self.bufs[5].ptr, self.bufs[6].ptr);
                let (mut n_chunks, mut k, mut chunk) = (chunk_read.len() as u32, k, CHUNK);
                let (mut bloom_mask, mut table_mask) = (self.bloom_mask, self.table_mask);
                let mut args: [*mut std::ffi::c_void; 13] = [
                    (&mut d_packed as *mut DevPtr).cast(),
                    (&mut d_starts as *mut DevPtr).cast(),
                    (&mut d_lens as *mut DevPtr).cast(),
                    (&mut d_chunk_read as *mut DevPtr).cast(),
                    (&mut d_chunk_first as *mut DevPtr).cast(),
                    (&mut n_chunks as *mut u32).cast(),
                    (&mut k as *mut u32).cast(),
                    (&mut chunk as *mut u32).cast(),
                    (&mut d_bloom as *mut DevPtr).cast(),
                    (&mut bloom_mask as *mut u64).cast(),
                    (&mut d_keys as *mut DevPtr).cast(),
                    (&mut table_mask as *mut u64).cast(),
                    (&mut d_hit as *mut DevPtr).cast(),
                ];
                let grid = n_chunks.div_ceil(BLOCK_THREADS);
                check(
                    (self.api.launch)(self.func, grid, 1, 1, BLOCK_THREADS, 1, 1, 0, std::ptr::null_mut(), args.as_mut_ptr(), std::ptr::null_mut()),
                    "kernel launch",
                )?;
                check((self.api.synchronize)(), "kernel")?;
                check((self.api.dtoh)(hit.as_mut_ptr().cast(), d_hit, n_reads), "copy from GPU")?;
            }
        }
        Ok(slot
            .into_iter()
            .map(|s| s.map(|r| r != u32::MAX && hit[r as usize] != 0))
            .collect())
    }
}

#[cfg(not(target_os = "linux"))]
pub struct Gpu;

#[cfg(not(target_os = "linux"))]
impl Gpu {
    pub fn open() -> Result<Self> {
        bail!("the GPU check is available on Linux only")
    }
    pub fn set_delta(&mut self, _: &BaitSet) -> Result<()> {
        unreachable!()
    }
    pub fn changed(&mut self, _: &[Entry]) -> Result<Vec<Option<bool>>> {
        unreachable!()
    }
}

fn as_bytes<T: Copy>(v: &[T]) -> &[u8] {
    // SAFETY: plain integer slices viewed as their bytes.
    unsafe { std::slice::from_raw_parts(v.as_ptr().cast(), std::mem::size_of_val(v)) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::CanonicalKmers;
    use crate::recruit::Rec;

    fn genome(len: usize, seed: u64) -> Vec<u8> {
        let mut x = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
        (0..len)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[(x >> 33) as usize & 3]
            })
            .collect()
    }

    fn entry(seq: Vec<u8>) -> Entry {
        let packed = seq.iter().all(|b| b"ACGT".contains(b)).then(|| {
            let mut p = vec![0u8; seq.len().div_ceil(4)];
            for (i, &b) in seq.iter().enumerate() {
                p[i / 4] |= (b"ACGT".iter().position(|&x| x == b).unwrap() as u8) << (2 * (i % 4));
            }
            p
        });
        Entry { rec: Rec { id: b"r".to_vec(), seq, qual: None }, packed }
    }

    /// Runs only where a GPU is (`OVASM_GPU=1 cargo test -- --ignored gpu`).
    #[test]
    #[ignore]
    fn gpu_hits_are_the_cpu_change_check() {
        let mut gpu = Gpu::open().expect("a CUDA GPU");
        for keys in [2_000usize, 600_000] {
            let mut bait = BaitSet::new(21);
            let t = bait.target("t").unwrap();
            bait.track_changes();
            let changed = genome(keys, 5);
            bait.add_sequence(t, &changed);
            let delta = bait.take_changes();
            gpu.set_delta(&delta).unwrap();
            let mut reads = vec![entry(Vec::new()), entry(b"ACGT".to_vec()), entry(genome(20, 3))];
            for i in 0..3000u64 {
                let len = 300 + (i as usize * 37) % 20_000;
                let mut s = genome(len, 100 + i);
                if i % 7 == 0 {
                    // one changed k-mer anywhere in the read, chunk edges included
                    let at = (i as usize * 13) % (len - 21);
                    let from = (i as usize * 101) % (keys - 21);
                    s[at..at + 21].copy_from_slice(&changed[from..from + 21]);
                }
                if i % 11 == 0 {
                    s[len / 2] = b'N';
                }
                reads.push(entry(s));
            }
            let got = gpu.changed(&reads).unwrap();
            let mut hits = 0;
            for (e, g) in reads.iter().zip(&got) {
                let cpu = CanonicalKmers::new(&e.rec.seq, 21).any(|(_, km)| delta.get(km) != 0);
                match g {
                    Some(g) => assert_eq!(*g, cpu),
                    None => assert!(e.packed.is_none()),
                }
                hits += usize::from(cpu);
            }
            assert!(hits > 300, "{hits} reads hold a changed k-mer");
        }
    }
}
