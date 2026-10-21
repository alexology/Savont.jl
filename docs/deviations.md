# Accepted deviations from upstream

Savont.jl reproduces savont 0.8.0 (commit `5a5a1ce`) exactly, except for the
upstream bugs listed here. Each deviation keeps an upstream-compatible mode so
the oracle comparison in `test/rust_oracle.jl` still checks exact parity.

## D1. First window of `minimizer_seeds_positions` (2026-10-21)

**Upstream behavior.** In `seeding::minimizer_seeds_positions`:

1. the first returned value is the last canonical k-mer of the first window as a
   raw 2-bit encoding, not its `mm_hash64` hash, paired with the position of the
   first window's minimum;
2. in the first window the forward k-mer is not masked to `2k` bits, so every
   k-mer after the first carries the preceding bases in its high bits, and its
   canonical choice and hash are wrong.

The main loop masks the k-mer and returns hashes, as does `fmh_seeds` in the same
file, so both are inconsistencies rather than design choices.

**Consequence.** Both callers use the values as content keys shared across
sequences: `merge.rs::fuzzy_merge_table` (export, `w = 28, k = 31`) merges an ASV
into a longer one, differing by at most 10 bp of terminal trimming, only if all
of the shorter ASV's minimizers occur in the longer one; and
`alignment.rs::remove_similar_seqs_kmers` (`w = 10, k = 21`) intersects candidate
sets over all minimizers. Because the first-window values depend on where the
sequence starts, a sequence trimmed at its 5' end never shares them with the
untrimmed one. On 50 simulated 1,300 bp amplicons, the trimmed copy's minimizer
set was contained in the original's in 0/50 cases with upstream behavior and in
50/50 with the correction whenever the 5' end was trimmed (3, 5 or 10 bp, with
`w = 28, k = 31` and `w = 10, k = 21`); with 3'-only trimming both gave 50/50.
Upstream therefore misses the 5'-trimmed pairs these functions are meant to find.

**Savont.jl.** `minimizer_seeds_positions` masks the first window and returns the
hash of the minimum as the first value, so every value is the hash of the
canonical k-mer at its reported position. Pipeline stages use this corrected
version. `upstream_quirks=true` reproduces upstream exactly and is used for the
oracle comparison. Expected effect when those stages are implemented: fuzzy
export merging and similar-consensus removal also catch 5'-trimmed sequences, so
results can contain fewer near-duplicate ASVs than upstream.
