# Differential validation plan

## Principle

Scientific equivalence is a measured property. The Rust binary built from commit `5a5a1ce9053a867abc09e991f61b8ec7fb4fc6ef` is the primary behavioral oracle. The preprint explains intent, while source and executable determine version-specific behavior.

Each run records Rust/Julia versions, package manifests, compiler flags, OS, CPU, thread count, command line, random seeds, input SHA256 hashes and reference commit. Inputs are immutable. Reports compare DNA sequences rather than arbitrary final ASV numbers.

## Oracle build

- Build upstream in release mode with its locked Cargo environment.
- Save the executable SHA256 and `cargo metadata`.
- Run Rust with 1, 2, 4 and available maximum threads to expose order sensitivity.
- Run the Julia reference mode with the same logical thread matrix.
- Preserve complete output directories and logs as CI artifacts for failed comparisons.

The Rust test data in `reference/savont/tests/data` is MIT-covered with the repository and contains two 1,000-read ONT FASTQ.gz files plus reference ASVs. Database downloads remain opt-in because of size and external licensing/availability.

## Instrumentation contract

The Rust oracle needs a small patch or companion harness that serializes stage data without changing algorithms:

1. filtered reads in stable input order, with length, estimated identity and sample index;
2. exact k-mer orientation counts, accepted SNPmers/blockmers and high-frequency set;
3. `TwinRead` sketches and marker vectors;
4. Stage-2 primary clusters as read-ID sets;
5. secondary clusters before and after every reclustering iteration;
6. ordered POA inputs and initial consensus;
7. quality-error map, pileup summaries, posterior calls and low-quality decisions;
8. merge equivalence classes and chimera parent decisions;
9. read-to-ASV candidate edges, likelihood/weight inputs and every EM iteration;
10. final per-sample depths and output tables.

Serialization order is explicit. Binary floating fields are accompanied by a decimal representation and units.

## Comparison rules

| Artifact | Comparator | Pass criterion |
|---|---|---|
| filtering, encodings, counts, markers | exact | byte/integer equality and identical ordered keys after canonical sort |
| clusters | set-based by original read ID | exact membership and exact number of clusters; representative separately exact |
| sequences | uppercase DNA, compare both orientation rules where applicable | exact sequence and boundaries for reference mode |
| pileup discrete observations | exact | exact counts per reference position and allele |
| probabilities/error rates | absolute + relative tolerance | initially `atol=1e-12`, `rtol=1e-10`; tighten or justify from operation order |
| EM | per-iteration vector after taxon/ASV key join | same convergence iteration where feasible; final `atol=1e-10`, integer depths exact |
| chimeras | exact query and parent sequences | exact flags, parent pair, match lengths and coverage decision |
| final ASVs | sequence-keyed set | exact set; headers checked after depth/order join |
| assignments | join by stable input file + record ordinal | exact ASV sequence or unassigned; identity tolerance `1e-8` percentage points |
| tables | schema and sequence-keyed rows | exact integer counts, taxonomy strings and column order |

No tolerance may convert a base, assignment, cluster, chimera flag or retained read from failure to success.

## Unit fixture matrix

### Input and qualities

- FASTA, FASTQ, gzip, multiline records and empty files;
- header whitespace and final `rc` token;
- absent qualities, minimum/maximum Phred, malformed lengths;
- reads exactly at length and 98% accuracy boundaries;
- pooled files with duplicate read names, proving ordinal identity is required.

### Encodings and markers

- all k-mers for small odd `k`, reverse complements and palindromes;
- ambiguous symbols at every window position;
- integer masks and wrapping hash vectors;
- central qualities 24/25/26;
- strand-count cells including zeros;
- binomial/Fisher decisions immediately below/at/above thresholds;
- high-frequency quantile ties and blockmer suffixes.

### Clustering

- identical reads;
- a single SNP at each amplicon region;
- two linked SNPs and conflicting unlinked errors;
- representative ties in accuracy and k-mer score;
- clusters at sizes 3/4 (HiFi) and 11/12 (default);
- no-SNPmer fractions around 10% and 75%;
- extreme cluster-depth ratios that trigger special merge branches;
- every tested input under multiple read permutations and thread counts.

### Consensus and posterior

- substitutions, insertions and deletions at internal and terminal sites;
- homopolymers of lengths 1–20 with and without HPC;
- quality-weight conflicts;
- fewer/equal/more than 75 POA reads and 250 pileup reads;
- length outliers affecting seed selection;
- band boundary versus `--no-band`;
- low-coverage end trimming and all low-coverage pileups;
- ambiguous-base masking and low-quality thresholds;
- pileup batch sizes 1, thread count, and full batch to prove scheduling invariance.

### Merge, chimera and EM

- true one-base ASVs at balanced and 1:10,000 abundance;
- errors one/three mismatches from a deep parent;
- exact two-parent chimeras with breakpoint sweeps;
- same-parent, highly similar-parent and insufficient-depth negative controls;
- unmapped reads, unique mappings and symmetric multimappings;
- all-zero and disconnected EM columns;
- pooled samples with shared/private ASVs and zero-count samples;
- final ordering ties and debug-to-final ID translation.

### Taxonomy and export

- miniature EMU, SILVA, GG2 and UNITE fixtures;
- species/genus identity threshold boundaries;
- taxonomic multimapping EM ties;
- exact SINTAX xorshift streams and bootstrap samples;
- SINTAX no-hit iterations and rank confidence denominator;
- exact, reverse-complement and terminal fuzzy ASV merging;
- conflicting lineages and duplicate sample labels.

## Synthetic biological scenarios

Generate reads from known haplotypes with a recorded simulator seed and truth table. Cross substitution rate, indel rate, homopolymer-context error, read accuracy, depth, abundance skew, amplicon length, strand balance and platform profile.

Required panels include:

- one to eight haplotypes separated by 1, 2, 5 and 20 variants;
- rare variants at 0.01%, 0.1%, 1% and 5%, with depths spanning the minimum cluster size;
- ONT R10.4 SUP-like and PacBio HiFi-like errors;
- 250 bp, full-length 16S, ~4 kb operon and variable ITS lengths;
- shared/private haplotypes across 2, 8 and 96 samples;
- explicit PCR chimeras at known breakpoints;
- protein-coding markers reserved for future extension tests, not reference acceptance.

For every panel report sensitivity, precision, sequence accuracy, abundance error, assignment accuracy, false chimera rate, runtime and peak RSS. These biological metrics supplement rather than replace exact Rust parity.

## Real-data matrix

1. Upstream ONT Zymo 1,000-read fixture, default parameters.
2. The second upstream Zymo fixture independently and pooled.
3. Combined upstream fixtures with thread determinism checks.
4. A licensed PacBio HiFi mock-community dataset.
5. A licensed ONT SUP environmental full-length 16S dataset.
6. An rRNA-operon or variable-length panel.

External datasets are added only with stable accession, checksum, license/terms and a small CI subset. Full datasets run in scheduled validation, not every pull request.

## Determinism

Repeat every small integration case at least three times for each thread count. Compare all stable output files byte-for-byte after normalizing log timestamps. Any nondeterminism is localized to the first divergent stage. Stable explicit sort keys are preferred over relying on dictionary or scheduler iteration order.

## Performance validation after parity

Benchmark reference and optimized modes against Rust after warm-up. Separate compilation latency, input decompression, computation and output. Record wall time, CPU time, peak RSS, Julia allocations, compressed bytes/s, reads/s and thread efficiency. Use identical storage and cold/warm cache policies.

Scale independently by read count, read length, number of ASVs and sample count. BenchmarkTools is appropriate for kernels; end-to-end runs use external wall-time/RSS collection in fresh processes. CPU profiling precedes any GPU proposal. A GPU experiment must report transfer time and end-to-end speed, not kernel time alone.

## Automated discrepancy report

Each differential run produces:

- `metadata.json` with provenance;
- one TSV/JSON comparison per stage;
- first divergent read/key/cluster plus surrounding evidence;
- aggregate exact/tolerance counts;
- final sequence-keyed ASV and abundance table;
- performance metrics;
- `report.md` linking failures to the relevant Rust and Julia functions.

CI fails on any unapproved discrete discrepancy. Accepted discrepancies live in a versioned registry containing rationale, impact, owner and expiry/review condition.
