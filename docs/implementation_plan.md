# Staged implementation plan

## Gate 0: architecture review

Review `architecture.md`, `algorithm_mapping.md`, this roadmap and `validation_plan.md`. Resolve the POA/minimap2 compatibility strategy and confirm whether the target is exact 0.8.0 behavior or a named later commit. No prototype result is promoted to reference status before this gate.

Deliverable: approved architecture record with the pinned upstream commit and explicit accepted discrepancies, initially none.

## Phase 2: package skeleton and observability

Create the final standalone `Savont.jl` repository beside the other Julia packages. Establish:

- standard `Project.toml`, `src/Savont.jl`, `test/runtests.jl`, `docs/`, `benchmark/`, `bin/` layout;
- upstream license attribution and a machine-readable `UPSTREAM.toml` containing repository, commit, version and license;
- immutable `SavontConfig` and preset tests;
- structured stage events with timings, counts, deterministic seeds and commit metadata;
- a reference/optimized mode boundary that prevents accidental algorithm substitution.

Acceptance: package loads on Julia 1.10 and current stable Julia on Linux and Windows; configuration/default snapshots match Rust CLI 0.8.0.

## Phase 3A: input, encoding and exact marker primitives

Implement streaming FASTA/FASTQ/FASTQ.gz, Cutadapt `rc` handling, quality pooling, packed DNA, reverse complements, two-bit rolling k-mers, `mm_hash64`, FracMinHash, split-k-mer masking and blockmers.

Use upstream functions as the oracle through small Rust fixture generators. Generate thousands of deterministic vectors, including ambiguous bases, palindromes, short reads, quality boundaries and `k=3:2:23`.

Acceptance: exact byte/bit/count equality for all primitive fixtures. Memory use must remain bounded during a multi-gigabyte synthetic stream.

## Phase 3B: count pipeline and SNPmer inference

Implement strand-separated first/second pass counting, optional Bloom prefiltering, exact count filters, high-frequency markers, binomial tests and Fisher strand-bias logic.

Start with a deterministic single-thread reference. Add multithreading only after equality. Parallel maps must be merged in explicit key order wherever order can propagate.

Acceptance: exact retained-read counts, k-mer counts, high-frequency sets and ordered `SnpmerInfo` records for synthetic fixtures and upstream Zymo data.

## Phase 3C: retained reads and hierarchical clustering

Implement `TwinRead` construction, length/accuracy filtering, SNPmer loading, low-polymorphism detection, primary FracMinHash clustering, secondary SNPmer/blockmer assignment and iterative consensus-marker reclustering.

Emit versioned diagnostic tables at every boundary. Do not implement clustering as generic edit distance or substitute DADA2 logic.

Acceptance: clusters compared as sets of original read IDs match Rust exactly at Stage 2, before reclustering and after convergence; iteration merge/reassignment counts also match.

## Phase 3D: consensus compatibility spike

Before writing the whole consensus stage, build a focused POA harness around `spoars` and candidate Julia approaches. Capture ordered input sequences, qualities, band settings, graph consensus and output from Rust.

Decision point:

1. if a Julia POA implementation matches all characterized behavior, use it;
2. otherwise isolate `spoars` behind a small compatibility interface for the reference mode while continuing a native Julia implementation as experimental;
3. document any platform/build consequences.

Then implement read selection/orientation, optional HPC, pileup generation, empirical quality-error estimation, posterior calls, homopolymer reconstruction, boundary trimming and low-quality filtering.

Acceptance: exact initial POA and final polished sequences for small/golden clusters; full Stage-4 diagnostics match or meet predeclared floating tolerances without changing base calls.

## Phase 3E: merge, chimera and quantification

Implement similar-consensus merging one decision rule at a time with fixtures for each branch. Implement the minimap2-backed chimera reference path, final candidate mapping, ASV EM, final sorting, assignments and pooled per-sample quantification.

Acceptance: exact final ASV sequence set, cluster membership, chimera flags and read assignments on upstream data; EM traces match under the tolerance policy; integer depths and per-sample feature tables match exactly.

## Phase 3F: taxonomy and export

Implement database registry/parsers, alignment classification EM, SINTAX's exact xorshift/resampling behavior, database downloads with checksum metadata, and exact/fuzzy QIIME export.

Acceptance: upstream integration tests are mirrored; local tiny EMU/SILVA/GG2/UNITE fixtures cover parsing without network; network database tests remain opt-in.

## Phase 4: full differential suite

Run Rust and Julia in controlled containers against the matrix in `validation_plan.md`. Produce machine-readable JSON/TSV discrepancies and a human Markdown report. Sequence-set comparison is canonicalized, while stage-level read/cluster IDs retain source identity.

The reference implementation exits Phase 4 only when all required gates pass. Any accepted deviation gets an issue, a documented biological consequence and a named compatibility status.

## Phase 5: profiling and Julia-specific optimization

Capture CPU profiles, allocations, peak RSS and I/O rates before changing algorithms. Candidate improvements, in likely order, are:

1. compact structure-of-arrays read/marker storage;
2. preallocated rolling encoders and thread-local counters;
3. measured Robin/Swiss dictionaries for integer keys;
4. deterministic chunked pipelines and disk-backed retained reads;
5. sparse read-to-ASV likelihood matrices and stable EM reductions;
6. parallel samples and clusters with explicit output ordering;
7. optional checkpoint schemas;
8. SIMD for fixed-width encoders and likelihood loops.

Every change runs the full differential suite and reports runtime, peak RSS, allocations and scaling. Optimizations live behind `OptimizedMode` until equality is established.

GPU experiments begin only when profiles identify a large regular kernel with enough arithmetic intensity and batch size. The likely experiments are batched encoding/count transforms or EM edge reductions. Hash-heavy discovery, variable-length POA and small-ASV workloads are poor first targets. KernelAbstractions may provide portability after a CPU/GPU cost model and transfer benchmark exist.

## Phase 6: extensions

Expose interfaces for error models, platform profiles, diagnostic confidence, replicate evidence, reading-frame validation and haplotype analyses. Extension output is explicitly namespaced and never changes reference-mode results.

## Dependency introduction policy

A dependency is added only when it supplies a maintained, tested primitive and reduces correctness risk or demonstrated engineering cost. Each proposed dependency gets:

- license and maintenance check;
- Linux/Windows and Julia 1.10 compatibility check;
- minimal fixture proving required semantics;
- performance/memory comparison for hot-path structures;
- a clear fallback or isolation boundary for native libraries.

The initial likely core set is FASTX, BioSequences, CodecZlib and Distributions. DataStructures is an evaluated performance dependency, not automatically selected. BenchmarkTools/Profile/Allocations support development only. GPU packages remain weak/extension dependencies.

## Immediate disposition of the existing prototype

The pre-architecture prototype contains algorithm substitutions (generic edit alignment/consensus, simplified clustering and chimera logic). It remains outside the standalone repository and the reference test target. Useful generic code may be promoted only after it passes the component gate associated with the corresponding Rust function.
