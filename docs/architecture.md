# SAVONT.jl architecture investigation

## Scope and reference baseline

This document covers Phase 1 only. It records the architecture to be reviewed before a reference implementation is written.

The inspected upstream checkout is:

- repository: <https://github.com/bluenote-1577/savont>
- branch at clone time: `main`
- commit: `5a5a1ce9053a867abc09e991f61b8ec7fb4fc6ef`
- commit date: 2026-10-10
- crate version: `0.8.0`
- upstream license declaration: `MIT OR Apache-2.0`; the repository currently contains an MIT `LICENSE`

The source checkout is preserved under `reference/savont`. The preprint is “Sensitive long-read amplicon sequence variant recovery with savont”, DOI `10.64898/2026.05.26.727271`. Its bioRxiv full text returned HTTP 403 during this investigation. Statements derived from the paper are therefore limited to searchable metadata and method summaries; implementation details below come from the pinned source code unless explicitly labelled otherwise.

Code drafted before this architecture review remains outside the standalone repository. It is an unvalidated prototype, is not a reference implementation, and must not be used as evidence of scientific equivalence.

### Investigation limits

The complete active Rust source, Cargo metadata, changelog, README and tests were inspected. The upstream Rust test suite was not executed because `cargo` is not installed on the current workstation; establishing a reproducible Rust oracle is consequently the first validation action after architecture review. The bioRxiv page for the preprint returned HTTP 403, so no claim below relies on inaccessible supplementary text. These are evidence gaps, not assumed passes.

## Active upstream boundary

`src/lib.rs` exports the active modules `kmer_comp`, `seeding`, `types`, `cli`, `constants`, `seq_parse`, `alignment`, `asv_cluster`, `chimera`, `classify`, `databases`, `download`, `merge`, `sintax`, `taxonomy`, and `utils`.

Several large files remain in the tree but are not compiled by version 0.8.0: `mapping.rs`, `map_processing.rs`, `polishing_mod.rs`, `polishing/`, `graph.rs`, and `unitig.rs`. They describe an earlier assembly/polishing architecture and contain useful tests and data structures, but copying them into the first reference implementation would reproduce dead code rather than SAVONT 0.8.0. This distinction must be retained in the module map and coverage reports.

## Reconstructed ASV workflow

### 0. CLI setup and presets

`main.rs::initialize_setup_cluster` validates files and the odd k-mer size, creates the output directory and Rayon pool, and applies presets:

- default full-length 16S: read length 1100–2000, minimum cluster size 12;
- `--hifi`: minimum cluster size 4;
- `--rrna-operon`: read length 3500–5000;
- default `k=17`, FracMinHash sampling factor `c=11`, minimum estimated read accuracy 98%, central-base quality 25;
- `--single-strand` changes strand-support requirements;
- `--use-hpc` enables homopolymer compression;
- `--pooled-samples` pools discovery but later quantifies each input separately.

Cutadapt reverse-complement annotations are recognized when the final read-header field is `rc`. Primer and adapter trimming is deliberately outside SAVONT.

### 1. Two-pass k-mer counting and polymorphic marker discovery

`seq_parse.rs::read_to_split_kmers` streams every input twice through bounded producer/consumer channels.

The first pass is optional when the requested Bloom-filter size is zero. Otherwise, strand-separated Bloom filters identify k-mers seen on both strands and reduce the exact counting set. Hash partitioning assigns k-mers to worker-local maps. The second pass exactly counts candidate split k-mers and experimental blockmers by orientation.

Candidates are retained when total multiplicity exceeds two and, unless single-strand mode is active, both orientations are observed. `kmer_comp.rs::get_snpmers_inplace_sort` then:

1. computes a high-frequency threshold from the upper tail of k-mer counts, with a minimum of 100;
2. sorts by the middle-base-masked split-k-mer context;
3. groups alternative central alleles at the same context;
4. orders alleles by abundance;
5. uses a binomial tail test to reject overly abundant secondary alleles compatible with sequencing error;
6. uses a two-sided Fisher exact test and odds ratio to reject strand-biased alternatives;
7. records accepted allele pairs as `SnpmerInfo`.

The exact threshold branches, zero-cell behavior, allele ordering and tie-breaking are compatibility requirements. The source currently calls `statrs` for the binomial CDF and `fishers_exact` for Fisher's exact test.

`kmer_comp.rs::twin_reads_from_snpmers` rereads the input, filters length and estimated quality, creates packed `TwinRead` objects, and annotates their sampled k-mers, SNPmers, blockmers, quality bins, input sample and orientation. Reads are sorted by estimated identity. If more than 75% have no SNPmers, low-polymorphism mode is enabled automatically; above 10% a warning is emitted.

### 2. Primary k-mer clustering

`asv_cluster.rs::cluster_reads_by_kmers` performs greedy clustering over reads sorted by estimated accuracy. It indexes representative FracMinHash k-mers and retrieves candidates by shared seeds. Candidate representatives are verified with prebuilt k-mer sets and the configured primary threshold. The stage preserves deterministic ordering explicitly where hash maps or parallel collections could otherwise change it.

After greedy assignment, small/compatible representatives are consolidated and clusters below the required size are removed. The output is written to `temp/kmer_clusters_stage2.tsv`.

This is not a generic edit-distance clusterer. Compatibility requires the upstream seed sampling, similarity calculation, representative order, candidate validation, merge rules and stable tie breaks.

### 3. Secondary SNPmer/blockmer clustering

`asv_cluster.rs::cluster_reads_by_snpmers` refines reads only within each primary k-mer group. It greedily assigns reads using compatible polymorphic marker alleles. The optional blockmer path changes marker extraction and comparison but retains the same hierarchy.

The stage then iterates consensus-marker reclustering:

- build a marker consensus for every subcluster;
- retain a marker allele when supported by at least `max(cluster_size / 6, 1)` reads;
- use median marker position and deterministic allele ties;
- merge clusters whose marker consensuses are concordant in both directions, with explicit exceptions for extreme depth ratios;
- reassign every read to the consensus with the fewest mismatches, then the most matches;
- remove clusters below `min_cluster_size`;
- repeat until both merges and reassignments converge or `max_iterations_recluster` (10) is reached.

Low-polymorphism mode skips this SNPmer split. Intermediate cluster tables are part of the differential-test surface.

### 4. POA consensus, pileup and quality filtering

`alignment.rs::align_and_consensus` runs per cluster in parallel. It selects at most 75 high-accuracy reads and a seed near the 90th percentile read length. A short-reference minimap2 index orients/maps selected reads to the seed. Sequences and quality values are optionally homopolymer-compressed and passed to native-Rust `spoars` overlap POA.

POA scoring is fixed at match `3`, mismatch `-8`, gap open `-6`, gap extend `-6`; the default alignment uses a band whose base is the maximum length deviation and whose fractional width is 0.1. `--no-band` uses exact unbanded overlap alignment. Input order, seed choice, quality encoding and the 75-read cap can change the consensus and must be matched.

All retained cluster reads are then mapped back to the POA consensus, up to `MAX_SEQS_CONSENSUS=250`, to build base/deletion/insertion pileups. Version 0.8.0 batches pileups under a memory budget. It builds only the deepest 10% first to estimate empirical errors by quality bin, drops those pileups, then regenerates and analyzes all clusters batch by batch.

`analyze_pileup_consensuses` trims ends below a cluster-specific coverage threshold and computes log posterior evidence for the reference and alternative at each site. Base qualities use the empirical quality-error map; indels use a fixed quality fallback. It records posterior diagnostics, optionally masks uncertain bases, estimates homopolymer lengths, and separates low-quality consensuses according to `lq_criteria`, including length, depth, ambiguous bases and posterior behavior.

The apparent separation between initial POA and pileup polishing is important: POA generates the candidate sequence; pileup analysis supplies empirical quality calibration, boundary correction, filtering and optional masking.

### 5. Similar-consensus merging

After HPC decompression, `alignment.rs::merge_similar_consensuses` considers both high- and low-quality consensuses. It uses minimap2 alignments, k-mer candidate indexes, depth relationships, mismatch counts, ambiguous-read evidence and deterministic equivalence classes to absorb sequences judged to be errors of deeper consensuses. Cluster membership and appended depth follow the surviving representative.

This stage is scientifically sensitive because an approximate replacement can either collapse true single-nucleotide ASVs or retain errors. It requires function-level characterization before implementation.

### 6. Chimera detection

`chimera.rs::detect_chimeras` first calculates pairwise consensus similarities. For every query it considers only potential parents with more than three times the query depth. It maps query to every potential parent using minimap2, extracts near-perfect left and right match lengths from CIGAR operations, handles reverse-complement coordinates, and evaluates parent-pair coverage and similarity restrictions. Combined parental segments must cover at least 95% of the query, subject to configurable allowable errors and minimum segment length (default `min_read_length / 10`). Detected queries are removed by `filter_chimeras`.

Pairwise mapping makes this a likely quadratic hot spot when ASV diversity is high.

### 7. Read-to-ASV mapping and EM depth refinement

`alignment.rs::refine_asv_depths_with_em` builds minimizer information for final candidates and maps retained reads to them. Version 0.8.0 precomputes candidate minimizer sets and uppercase sequences. It uses identity thresholds `[0.99, 0.9975, 1.0]` with equal coverage weights, includes low-polymorphism minimap2 fallback behavior, and records raw mappings in `temp/read_to_asv_mappings.tsv`.

EM distributes ambiguous reads among compatible ASVs. Candidate weights and alignment evidence, convergence criteria, rounding/depth conversion and all-zero behavior must be extracted into golden unit cases. Survivors are sorted by refined depth, then sequence, then original ID. A debug-ID-to-final-index map is retained so final assignments remain correct after sorting.

For pooled runs, `compute_per_sample_depths` repeats quantification by `sample_idx` against the same discovered ASV set. The sorted EM FASTA is rewritten before this step; upstream 0.6.3 fixed an index/order bug here, making sort order an explicit invariant.

### 8. Outputs

The ASV command writes:

- `final_asvs.fasta`, with `final_asv_N_depth_M` identifiers;
- `feature-table.tsv`, an ASV-by-sample QIIME2-compatible count table;
- `final_assignments.tsv`, including unassigned retained reads, alignment identity and estimated read identity;
- `temp/final_clusters.tsv` and other stage-specific diagnostics.

Output comparison must join by DNA sequence, not `N`, because IDs are derived from a final ordering.

## Other command workflows

### Alignment classification

`classify.rs` loads ASVs and depths, builds or reuses a minimap2 database index, keeps best mappings, and runs a separate EM over taxonomic targets to resolve multi-mapping. Its convergence threshold is `0.1 / total_reads`. The target with greatest EM abundance determines the principal assignment, while sequence-identity thresholds control species/genus labels. Pooled feature-table counts are carried into taxonomic abundance outputs.

### SINTAX classification

`sintax.rs` extracts canonical 12-mers. For every ASV and iteration it deterministically seeds a custom xorshift RNG and samples 32 k-mers with replacement. A reverse query map lets reference sequences be streamed in parallel; reference k-mers are deduplicated, and the best-scoring taxonomy per bootstrap is retained. Rank confidence is the number of votes divided by total iterations, including iterations without a classification in the denominator. Species output is intentionally suppressed: this command is genus-level maximum.

### Export

`merge.rs` canonicalizes exact sequence hashes with reverse complements, combines feature tables and taxonomy, and optionally merges sequences differing only by terminal trimming of at most 10 bp. It writes QIIME2 feature, representative sequence and optional taxonomy tables. Deterministic hashes, representative selection, sample relabeling, fuzzy merge direction and lineage conflict behavior are test targets.

## Proposed Julia package architecture

The reference implementation should use small modules with an explicit data flow:

```text
Savont
├── Configuration        # immutable CLI/reference configuration and presets
├── IO                   # streaming FASTA/FASTQ(.gz), output schemas, checkpoints
├── Sequences            # packed DNA, qualities, reverse complement, HPC
├── Kmers                # exact encodings, hashing, split/block/SNP markers
├── ReadStore            # compact retained-read storage and sample provenance
├── Clustering
│   ├── Primary          # FracMinHash representative clustering
│   └── Polymorphic      # hierarchical SNPmer/blockmer clustering and reclustering
├── Consensus
│   ├── POA              # reference-compatible POA interface
│   ├── Pileup           # read mapping and bounded-memory pileups
│   └── Posterior        # quality-error estimation and consensus decisions
├── Refinement           # similar-consensus merging and chimera filtering
├── Quantification       # read mappings, ASV EM and per-sample depths
├── Taxonomy             # database registry, alignment classification and SINTAX
├── Export               # QIIME2 and multi-run merging
├── Diagnostics          # structured stage events and intermediate snapshots
└── CLI                  # thin argument-to-configuration adapter
```

Reference semantics should be selected by a `ReferenceMode` dispatch marker. Later optimized implementations may implement the same narrow interfaces under `OptimizedMode`, but may not silently replace the reference path.

Key concrete types should include packed `ReadRecord`, `TwinRead`, `SnpmerInfo`, `ClusterSet`, `ConsensusSequence`, `Pileup`, `ReadMapping`, `EMState`, `TaxonomyEntry`, and stage-specific diagnostic records. IDs should be integer indices into immutable stores rather than object references. Storage layout should be chosen only after measuring the Rust baseline and Julia prototypes.

## Expected computational hot spots

Based on source structure and the 0.8.0 changelog, the likely dominant components are:

1. two complete FASTQ passes and k-mer extraction/counting;
2. primary candidate retrieval and representative verification;
3. POA for each retained cluster;
4. read-to-consensus minimap2 mapping and large pileups;
5. similar-consensus merging and pairwise chimera mappings;
6. final read-to-ASV mapping and EM;
7. reference-database indexing/classification for large databases.

These rankings are hypotheses until profiles are collected. The changelog reports upstream improvements of about 3.3× for clustering, 1.6× for EM refinement and up to 2.9× lower Stage-4 peak memory, which identifies the same pressure points.

## Julia opportunities and constraints

Julia can provide concrete packed types, generated/specialized k-mer code, thread-local counters, composable iterators, multiple-dispatch backends, memory-mapped intermediates, sparse EM matrices and direct numerical kernels without an FFI boundary. Existing local packages contribute tested conventions: DADA2.jl already uses `BioSequences`, `CodecZlib` and `Distributions`; ProfileAlignments.jl contains Julia pair/profile alignment code useful for experiments and tests.

ProfileAlignments must not become an implicit reference dependency: its algorithms differ from `spoars` overlap POA and minimap2. Reuse is acceptable only behind a characterized interface after differential validation.

GPU work is deferred. K-mer extraction and sparse reductions might be suitable at high batch sizes, but irregular hash tables, POA graphs, variable-length alignments and PCIe transfers can erase benefits. GPU code enters only after reference validation and CPU profiling demonstrates a realistic end-to-end gain.

## Open questions for Phase 2

- Whether `spoars` behavior can be reproduced exactly in Julia or must be isolated behind a temporary native compatibility backend.
- Which minimap2 details influence results beyond documented presets and bucket sizing.
- Exact semantics of every branch in similar-consensus merging and ASV EM.
- Whether Bloom filtering is active in common default runs (`bloom_filter_size=0`) and how auto-sizing is intended externally.
- Which intermediate files are stable interfaces versus debugging artifacts.
- How the inaccessible preprint describes thresholds or assumptions not obvious from code.
