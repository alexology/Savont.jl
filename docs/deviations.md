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
results can contain fewer near-duplicate ASVs than upstream. Measured on the export
test runs (`test/data/taxonomy/export_*`, 11 distinct sequences): upstream merges
only the copy trimmed by 5 bases at the 3' end (10 ASVs), the correction also the
copies trimmed by 4 and 2 bases at the 5' end (8 ASVs).

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

## D3. Low-polymorphism primary clusters (2026-10-21, methodological change)

This is not an upstream coding error but a change of method requested by the
user.

**Upstream behavior.** In stage 3, a read without SNPmers matches no
representative in the greedy SNPmer pass, starts a singleton cluster and is
removed by `min_cluster_size`. Upstream avoids this only through its
low-polymorphism mode, enabled by the user or automatically when more than 75%
of all reads lack SNPmers, which then skips stage 3 for every cluster. An
amplicon without SNPmers in a sample whose other amplicons are polymorphic is
therefore discarded entirely: on the test data, an unrelated amplicon of 45
reads disappears. In metabarcoding this can affect, for example, a species
represented by a single haplotype and without close relatives in the sample.

**Savont.jl.** The 75% rule is applied to each primary cluster: a primary
cluster in which more than 75% of the reads have no SNPmer skips stage 3 and is
kept whole (if it has at least `min_cluster_size` reads), as the whole data set
is in low-polymorphism mode. Keeping only the SNPmer-less reads together inside
polymorphic primary clusters was not chosen: there they are mostly reads whose
SNP sites had low quality, so they would form a cluster mixing several
haplotypes and a mixed consensus. `SavontConfig(upstream_quirks=true)`
restores the upstream behavior. On the test data the 45-read amplicon is kept
and the other clusters are unchanged.

## D4. Extra base at the end of the initial consensus (2026-10-21)

**Upstream behavior.** `utils::homopolymer_compress(seq, do_hpc)` with
`do_hpc = false` copies the sequence base by base, but the push of the final run
after the loop is outside the `do_hpc` branch, so it appends `seq[0]` (and a
length of 1) at the end. The only caller is `alignment.rs` (stage 4a), which
passes the spoars consensus through it with `do_hpc = use_hpc`; without
homopolymer compression (the default) every initial consensus therefore ends
with a copy of its first base. Returning the input unchanged is plainly the
intent of `do_hpc = false`, and the compressed branch handles the final run
correctly.

**Consequence.** On the test data, all 38 initial consensuses of the 11 oracle
configurations are one base longer than the POA consensus, and the extra base is
the first base. In stage 4b the extra base is almost never covered by reads, so
it falls in the low-coverage end that is masked with `N` and removed on
decompression: polishing the same clusters with and without it gave identical
decompressed consensuses and low-quality calls in 83 of 83 cases (test data and
20,000 real ONT reads of a 2.8 kb amplicon). The correction is therefore
cosmetic for the pipeline; it matters for direct users of `initial_consensus`.

**Savont.jl.** `homopolymer_compress(seq, false)` returns the sequence
unchanged. `upstream_quirks=true` appends the first base as upstream and is used
for the oracle comparison.

## D5. k-mer deduplication of consensuses never removes anything (2026-10-21)

**Upstream behavior.** `alignment.rs::remove_similar_seqs_kmers` (start of stage 5)
is meant to drop a consensus whose minimizers all occur in a consensus with more
than twice its depth ("Remove subsetted reads at lower depth"). It keeps a
consensus by moving it out of the vector with `std::mem::take`, which leaves a
default value with depth 0 in its place, and later comparisons read the depth
from that vector. The consensuses arrive sorted by decreasing depth, so every
deeper candidate has already been kept, and zeroed, when a consensus is
examined, while the later ones are not deeper: no consensus is ever removed. The
only remaining effect is that consensuses shorter than 100 bases are dropped.

**Consequence.** On 100,000 real ONT reads of a 2.8 kb amplicon, the corrected
deduplication removes 8 of 181 polished consensuses (upstream: none). The
alignment-based merge that follows would have merged all 8 anyway: the final
stage-5 consensuses are the same 76 sequences in both modes, but the reads of the
removed consensuses (1,989 of 61,645) are no longer added to the clusters they
would have merged into. On the 20,000-read subset and on the simulated data the
results are identical, except for the 2%-error simulated reads, where a 30-read
consensus is removed instead of merged. Stage 7 re-assigns reads to the final
ASVs, so the effect on reported depths is expected to vanish there (to be
verified when stage 7 is implemented).

**Savont.jl.** `remove_similar_consensuses` compares the original depths.
`upstream_quirks=true` compares against zero for consensuses already kept, as
upstream.

## D6. Chimera detection never finds the pairwise similarities (2026-10-21)

**Upstream behavior.** `chimera.rs::calculate_pairwise_similarities` maps each
consensus j to each earlier consensus i (i < j) and stores the identity under the
key `(j, i)`. Every lookup in `detect_chimeras` uses `(min, max)`, i.e. `(i, j)`,
so none is ever found and the defaults apply: 1.0 between a query and a parent,
0.0 between the two parents. Hence: the chimera score is always 0; the
parent-similarity condition (parents below 97% similar, or below 99.5% if both
are more than 10 times deeper) always holds and the coverage threshold
min(0.9 max(s, 0.7), 0.8) is always 0.63; the second rule (a single parent
covering all but the last bases with more than 20 / log2(depth ratio) estimated
mismatches) never fires, because the estimated mismatches are always 0.

**Consequence.** These rules were written but never active upstream, so its
thresholds were tuned without them; the correction is therefore a change of
results whose quality is not established. On the simulated chimera data
(`test/data/chimeras.fastq.gz`), the correction keeps a true chimera whose
parents are 98.5% similar and less than 10 times deeper (upstream removes it), and
removes a true variant with 5 substitutions in its last 100 bases (single-parent
rule; upstream keeps it). On 20,000 real ONT reads (42 consensuses after stage 5,
22 called chimeric upstream), the correction keeps one consensus of depth 20
whose two "parents" are 99.76% similar and whose left and right matches together
cover 159% of it, i.e. a close variant of both rather than a chimera. On
100,000 real reads (76 consensuses after stage 5, everything else as upstream),
upstream calls 51 chimeras and reports 25 ASVs, the correction 48 and 28: the
three extra ASVs (1,094, 193 and 138 reads) differ by 2, 10 and 15 bases from
ASVs that upstream reports deeper by about the same number of reads. This is the
only deviation that changed the final ASVs on these data.

**Savont.jl.** `detect_chimeras` looks up the stored pairs. `upstream_quirks=true`
looks them up as upstream. The effect is mixed; the user reviewed it and chose
the correction as the default (2026-10-25). Upstream still has the bug at its
latest commit (`5a5a1ce`, checked 2026-10-25): the pairs are stored as `(j, i)`
with `j > i` in `calculate_pairwise_similarities` and looked up as `(min, max)`.

## D7. Stage 5 panics when no consensus is 100 bases long (2026-10-21)

**Upstream behavior.** `alignment.rs::merge_similar_consensuses` returns early
only when its input is empty. `remove_similar_seqs_kmers` then drops the
consensuses shorter than 100 bases; if none is left, the empty FASTA it writes
makes the minimap2 crate return `Err("Index File is empty")`
(`Aligner::with_index`), and `.expect("Failed to create aligner")` panics.

**Consequence.** Amplicons shorter than 100 bases (after trimming) stop the
whole run instead of giving no ASV.

**Savont.jl.** `merge_similar_consensuses` returns an empty vector.
`upstream_quirks=true` throws an error, as upstream stops.

## D8. `--skip-chimera-detection` ends the run (2026-10-21)

**Upstream behavior.** In `main.rs::run_cluster`, the option is handled with
`if args.skip_chimera_detection { log::info!("Skipping chimera detection as per
user request."); return; }` right after stage 5, which returns from the whole
command: stage 7 (EM quantification) and the output files (`final_asvs.fasta`,
`feature-table.tsv`, `final_assignments.tsv`) are never produced, although the
message says only chimera detection is skipped.

**Consequence.** With the option, upstream leaves only intermediate files.

**Savont.jl.** `savont_asv` skips stage 6 only and completes the run.
`upstream_quirks=true` stops after stage 5 like upstream (`asvs = nothing`, no
output files).

## D9. Reads with tab-separated header tags are reported unassigned (2026-10-21)

**Upstream behavior.** A read's id is its whole FASTQ header. Stage 7 writes one
row per assignment to `temp/read_to_asv_mappings.tsv` (`id`, `debug_id:N`, ...),
and `alignment.rs::write_final_assignments` reads it back with
`line.splitn(5, '\t')`. When the header contains tabs, as in the FASTQ written by
ONT's dorado (`qs:f:...`, `RG:Z:...` tags after the read name), the second field
is a tag instead of `debug_id:N`, the row is ignored, and the read is written as
`unassigned`; its header, tabs included, also shifts the columns of
`final_assignments.tsv`.

**Consequence.** On real dorado reads (20,000 reads of a 2.8 kb amplicon), every
read of `final_assignments.tsv` is unassigned and the file has a variable number
of columns, while the ASVs and depths are unaffected.

**Savont.jl.** Assignments are passed in memory, and the `read_id` column holds
the header up to its first tab (the read name). Headers without tabs give the
same file as upstream. `upstream_quirks=true` reproduces upstream's file.

## D10. A gzip-compressed SILVA taxonomy map cannot be read (2026-10-22)

**Upstream behavior.** `taxonomy.rs::Database::load_silva` looks for the taxonomy
map as `taxmap_*.txt` or `taxmap_*.txt.gz`, but `load_silva_taxonomy` reads the
file as plain text (`BufReader::new(File::open(path))`), while the FASTA and the
GreenGenes2 headers next to it are decompressed when their name ends in `.gz`.

**Consequence.** A SILVA directory holding the map as distributed
(`taxmap_slv_ssu_ref_nr_138.2.txt.gz`) stops `classify` and `sintax` with
`stream did not contain valid UTF-8` (checked with the oracle on
`test/data/taxonomy/silva-138.2` with the map compressed). Directories made by
`savont download` are not affected, because it decompresses the map.

**Savont.jl.** `load_database` decompresses the map; the resulting database is
the same as with the plain file. `upstream_quirks=true` fails as upstream.

## D11. An even sampling rate separates the two strands (2026-10-22)

**Upstream behavior.** `get_twin_read_syncmer` keeps a k-mer as an open syncmer
when the smallest of its `c` s-mer hashes is at index `(c - 1) / 2` (integer
division). For odd `c` this is the central s-mer, the same on both strands; for
even `c` there are two central s-mers, and the reverse complement swaps them, so
a k-mer selected on one strand is never selected on the other. `-c` accepts any
value up to `k`.

**Consequence.** On 1,300 random bases and their reverse complement, `c = 11`
selects the same syncmers on both strands, `c = 10` none in common (also
`c = 2`: 638 and 646 syncmers, none shared; review SV1), so reads of the two
strands of one amplicon fall in different primary clusters. Defaults and presets
use `c = 11`.

**Savont.jl.** `validate` rejects even `sampling_rate`. `upstream_quirks=true`
accepts it, as upstream.

## D12. EM depths are rounded one by one (2026-10-22)

**Upstream behavior.** Stage 7 (`refine_asv_depths_with_em`, and per sample in
`compute_per_sample_depths`) sets each depth to `round(abundance × assigned)`
independently.

**Consequence.** The depths need not add up to the assigned reads: one read
ambiguous between two ASVs gives depths 1 and 1, between three ASVs 0, 0 and 0
(the read is lost and the ASVs may be dropped as depth 0), 31 reads between two
ASVs 16 and 16 (review SV2). On 20,000 real ONT reads, one ASV of 21 changes
from 302 to 303 reads; ASVs, sequences and read assignments are unchanged, and
the total equals the 13,116 assigned reads.

**Savont.jl.** Largest-remainder rounding: integer parts, plus one read for the
largest fractional parts (earlier ASV on ties) until the assigned reads are
reached; every depth differs by at most one from upstream's.
`upstream_quirks=true` rounds as upstream.

## D13. Reads sharing an id share their assignment (2026-10-22)

**Upstream behavior.** `write_final_assignments` reads `read_to_asv_mappings.tsv`
back into a map keyed by read id; for an id seen twice (the same read name in two
pooled input files, for example), the first read's ASV and NM are written for
both.

**Consequence.** Wrong `asv` and `alignment_identity` columns for the later
reads with a repeated name (review SV3); depths are unaffected.

**Savont.jl.** Each read is written with its own first mapping row.
`upstream_quirks=true` matches by id as upstream.

## D14. Results of an earlier run survive an empty rerun (2026-10-22)

**Upstream behavior.** `export` writes `merged_taxon_counts.tsv` only when there
is at least one ASV and does not remove an older one; `sintax` on an empty
`final_asvs.fasta` stops with an error and writes nothing.

**Consequence.** Rerunning into the same directory can leave a taxon table, or
sintax tables, of an earlier run next to the new files (review SV4).

**Savont.jl.** `savont_export` removes an older `merged_taxon_counts.tsv` when
there is no ASV; `savont_sintax` writes both tables with their header only. Only
file sets of a rerun differ; no mode reproduces the stale files.

## D15. Classify divides by a total depth of 0 (2026-10-22)

**Upstream behavior.** `classify` takes depths from `feature-table.tsv`, with 0
for ASVs whose name is not in it, and divides by their total without a check.

**Consequence.** With a table of zeros, or one whose ids do not match the FASTA
(for instance a table from another run), every abundance is `NaN` in the output
tables (review SV5).

**Savont.jl.** `savont_classify` warns about ASVs missing from the table and
throws an `ArgumentError` when the total depth is 0. `upstream_quirks=true`
writes the `NaN` abundances as upstream.

## D16. Pooled export puts unmatched ASVs in the first sample (2026-10-25)

**Upstream behavior.** `merge.rs::export` takes an ASV's depths from the row of
`feature-table.tsv` with the same name; for an ASV without a row it falls back to
the total of the depths in its FASTA header, placed in the run's first sample
column (`merge.rs:454–457`).

**Consequence.** For a pooled run whose feature table does not match its FASTA
(edited or from another run), the reads of every sample are attributed to the
first one without any message: header `final_asv_0_depth_3-7` with samples S1 and
S2 gives `[10, 0]` instead of `[3, 7]` (review SV7, reproduced with the oracle).

**Savont.jl.** `savont_export` warns about ASVs missing from a run's table and
uses the per-sample depths of a pooled header (one per sample column); a header
that does not give them is an error. Runs without a feature table, and
single-sample tables, are unchanged. `upstream_quirks=true` falls back as
upstream (the warning is still given).

## Upstream behavior kept on purpose

- **Reads without SNPmers inside polymorphic primary clusters** are still
  discarded in stage 3, as upstream (see D3).
- **Banded partial order alignment depends on the vector width.** spoars, used by
  upstream for the initial consensus, computes banded rows in segments of 16 or 8
  columns (AVX2, 16- or 32-bit scores) or 8 or 4 (SSE4.1, NEON), so the same
  reads can give a different consensus on different CPUs when the optimal
  alignment leaves the band. `poa_align` reproduces both (`isa = :avx2`, the
  default, or `:sse41`).
- **Reclustering stops when a round makes no merge**, although the upstream
  comment mentions reassignments too.
- **Stage 4b edge cases** (`analyze_pileup!`): when no pileup position reaches
  the minimum coverage, nothing is trimmed (the warning branch for this case is
  unreachable upstream); `right_end - 100` is an unsigned subtraction that wraps
  in release builds, so consensuses shorter than 100 positions get no right end
  region; inserted bases contribute only through the base before them (their own
  terms are taken with `.take(0)`); the left end mask stops before the
  low-confidence position that defines it, the right one includes it. None is a
  demonstrable error, and all affect only short or poorly covered consensuses.
- **Merging of low-quality consensuses** (stage 5): upstream maps them to the
  other consensuses and adds their depth to an `appended_depth` field, which is
  reset when the merged consensuses are rebuilt a few lines later, so the step
  has no effect. `merge_similar_consensuses` omits it; the result is the same.

- **High-frequency filtering of syncmers and SNPmers** (`kmer_comp::twin_reads_from_snpmers`)
  removes entries from `minimizer_positions` and `snpmer_positions` only; the
  lines that would filter `minimizer_kmers` and `snpmer_kmers` are commented out
  upstream, so primary clustering and the LSH signatures still use every
  syncmer. This looks deliberate and is reproduced (see `TwinRead`).
- **Taxonomy subcommands** (classify, sintax, export) reproduce upstream's
  choices, including: classify keeps the alignments with the NM of the primary
  one (not the smallest NM), measures identity over the aligned query span, and
  among taxa of equal EM abundance chooses the last alignment while
  `asv_mappings.tsv` lists the first one first (so `export`, which takes the
  first row, can report the other taxon); sintax counts a k-mer drawn twice
  twice and never assigns species; the pooled genus table merges lineages that
  differ only in the superkingdom; reference names are cut at the first
  whitespace for alignment (GreenGenes2 headers with a space are then never
  found by classify) but not for sintax.
- **Upstream nondeterminism made deterministic.** Upstream writes rows of equal
  abundance of the single-sample abundance tables in hash-map order, breaks ties
  between equally frequent sintax names, and between equally deep candidates of
  the export fuzzy merge, in hash-set order, and resolves ties between references
  of equal sintax score by thread timing. Savont.jl writes tied rows in order of
  first appearance, takes the first name seen, the smallest feature ID and the
  earliest reference (upstream with one thread). The oracle comparisons run
  upstream with one thread and sort tied rows; the synthetic cases have no ties
  of names or merge candidates (checked by repeated oracle runs).
