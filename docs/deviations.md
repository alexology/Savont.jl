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

## D2. Syncmers at the start of each read (2026-10-21)

**Upstream behavior.** `seeding::get_twin_read_syncmer` selects open syncmers
from a rolling window of `k - s + 1` canonical s-mer hashes (`s = k - c + 1`, 7
with the defaults `k = 17, c = 11`). Its initialization loop adds only the first
`s - 1` bases to the rolling s-mer, and the main loop then continues from base
`k - 1`, skipping bases `s - 1` to `k - 2`. Consequently:

1. no hash is computed for the s-mers inside the first k-mer, so the k-mers
   starting at positions `0` to `k - s - 1` (0–9 with the defaults) are never
   examined and can never be selected;
2. the next `s - 1` hashes are computed from non-contiguous bases, so the
   k-mers starting at positions `k - s` to `k - 2` are judged on wrong hashes.

The syncmer property of a k-mer depends only on its sequence; these choices
depend on its position in the read, so they are an initialization error rather
than a design choice.

**Consequence.** Read starts vary because of primer trimming and indels. On 100
simulated 1,300 bp reads (`k = 17, c = 11`), removing the first 3 bases changed
the syncmers of the remaining sequence beyond the 3-base shift in 76/100 reads
with upstream behavior and in 0/100 with the correction; the earliest syncmer was
at position 10 or later in every read with upstream behavior (median 12) and at
median position 2 with the correction. Syncmers feed primary clustering (shared
minimizers and the LSH bucket signatures), so reads of the same ASV that start at
different positions share fewer markers than they should. On the test data the
final clusters are identical in both modes in 7 of 11 configurations; the others
(stricter primary threshold, `c = 7`, simulated reads with 2% substitutions)
differ by a few reads per cluster or by one extra small cluster, with no
systematic gain or loss in cluster purity against the simulated haplotypes.

**Savont.jl.** The rolling s-mer is filled with every base before the first
k-mer, so every k-mer of the read is examined. `SavontConfig(upstream_quirks=true)`
reproduces upstream exactly and is used for the oracle comparison.

## Upstream behavior kept on purpose

- **Reads without SNPmers are discarded in stage 3** unless low-polymorphism mode
  is on (set by the user or automatically when more than 75% of the reads lack
  SNPmers): in the greedy SNPmer pass such a read matches no representative, starts
  a singleton cluster and is removed by `min_cluster_size`. On the test data an
  unrelated amplicon (45 reads, no SNPmers) is lost entirely. Upstream warns when
  more than 10% of the reads lack SNPmers and its low-polymorphism comment states
  that SNPmer-less reads are otherwise discarded, so this is a known design
  limitation rather than a coding error.
- **Reclustering stops when a round makes no merge**, although the upstream
  comment mentions reassignments too.

- **High-frequency filtering of syncmers and SNPmers** (`kmer_comp::twin_reads_from_snpmers`)
  removes entries from `minimizer_positions` and `snpmer_positions` only; the
  lines that would filter `minimizer_kmers` and `snpmer_kmers` are commented out
  upstream, so primary clustering and the LSH signatures still use every
  syncmer. This looks deliberate and is reproduced (see `TwinRead`).
