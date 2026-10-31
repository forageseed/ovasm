// The replay passes' check for changed k-mers (`gpu.rs`), one thread per chunk of a read's
// k-mer starts. The lookup is `bait.rs`'s `LookupIndex` bit for bit: the same hash, filter bits,
// table slot and linear probing, so a read is hit here exactly when `BaitSet::any_key` says so.
// Compiled to PTX by `scripts/build_gpu_ptx.py`; the PTX is committed as `gpu_delta.ptx`.

typedef unsigned long long u64;
typedef unsigned int u32;
typedef unsigned char u8;

__device__ __forceinline__ u64 mix(u64 key) {
    u64 h = (key ^ (key >> 33)) * 0xFF51AFD7ED558CCDull;
    return h ^ (h >> 33);
}

__device__ __forceinline__ bool bit(const u64* block, u32 b) {
    return (block[b >> 6] >> (b & 63)) & 1;
}

// Reads are 2-bit packed (base i at bits 2(i % 4) of byte i / 4) at `packed + starts[r]`.
// Chunk c covers the k-mers starting at chunk_first[c] .. chunk_first[c] + chunk of read
// chunk_read[c]; any hit sets hit[r].
extern "C" __global__ void any_delta(
    const u8* packed, const u64* starts, const u32* lens,
    const u32* chunk_read, const u32* chunk_first, u32 n_chunks, u32 k, u32 chunk,
    const u64* bloom, u64 bloom_mask, const u64* keys, u64 table_mask, u8* hit)
{
    u32 c = blockIdx.x * blockDim.x + threadIdx.x;
    if (c >= n_chunks) return;
    u32 r = chunk_read[c];
    const u8* p = packed + starts[r];
    u32 len = lens[r];
    u32 first = chunk_first[c];
    u32 end = min(first + chunk, len - k + 1);
    u64 mask = (1ull << (2 * k)) - 1;
    u32 shift = 2 * (k - 1);
    u64 fwd = 0, rev = 0;
    for (u32 i = first; i + 1 < first + k; i++) {
        u64 code = (p[i >> 2] >> (2 * (i & 3))) & 3;
        fwd = ((fwd << 2) | code) & mask;
        rev = (rev >> 2) | ((3 - code) << shift);
    }
    for (u32 s = first; s < end; s++) {
        u32 i = s + k - 1;
        u64 code = (p[i >> 2] >> (2 * (i & 3))) & 3;
        fwd = ((fwd << 2) | code) & mask;
        rev = (rev >> 2) | ((3 - code) << shift);
        u64 key = fwd < rev ? fwd : rev;
        u64 h = mix(key);
        const u64* block = bloom + ((h & bloom_mask) << 3);
        if (bit(block, (h >> 37) & 511) && bit(block, (h >> 46) & 511) && bit(block, (h >> 55) & 511)) {
            u64 slot = ((h << 32) | (h >> 32)) & table_mask;
            for (;;) {
                u64 stored = keys[slot];
                if (stored == key) { hit[r] = 1; return; }
                if (stored == ~0ull) break;
                slot = (slot + 1) & table_mask;
            }
        }
    }
}
