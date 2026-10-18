# Rust-to-Julia algorithm mapping

This map targets upstream commit `5a5a1ce9053a867abc09e991f61b8ec7fb4fc6ef` (savont 0.8.0). Status values are `reference-required`, `support`, `deferred`, and `inactive-upstream`.

| Rust source / responsibility | Proposed Julia module | Status | Compatibility notes |
|---|---|---:|---|
| `main.rs` stage orchestration | `Pipeline` | reference-required | Preserve stage order, sorting, file boundaries and pooled stage 7b |
| `cli.rs`, presets and defaults | `Configuration`, `CLI` | reference-required | CLI is a thin adapter; immutable config is the tested API |
| `constants.rs` | `Configuration.Constants` | reference-required | Values versioned with upstream commit; avoid unexplained magic literals |
| `seq_parse.rs` streaming, two-pass count, Bloom filters | `IO.Readers`, `Kmers.Counting` | reference-required | Preserve `rc` header handling, orientation counts, thresholds and partition determinism |
| `seeding.rs` two-bit encodings, hashes, FracMinHash, split/block markers | `Sequences.Encoding`, `Kmers.Seeding` | reference-required | Bit-level golden vectors and overflow semantics are mandatory |
| `kmer_comp.rs` SNPmer statistics, retained-read construction, HPC | `Kmers.Polymorphism`, `ReadStore`, `Sequences.Homopolymers` | reference-required | Match binomial/Fisher tails, tie order, quality bins and read sort |
| active portions of `types.rs` | `CoreTypes` split among owning modules | support | Do not port historical unused graph types merely because they remain public in Rust |
| `asv_cluster.rs` primary k-mer clustering | `Clustering.Primary` | reference-required | Candidate index, representative ordering and validation are algorithmic behavior |
| `asv_cluster.rs` SNPmer/blockmer refinement | `Clustering.Polymorphic` | reference-required | Preserve hierarchy, depth exceptions, marker consensus and convergence |
| `alignment.rs::generate_consensus_poa` via `spoars` | `Consensus.POA` | reference-required | Highest fidelity risk; compare graph/consensus behavior on exhaustive small cases |
| `alignment.rs` minimap2 orientation and pileup generation | `Consensus.Pileup` | reference-required | Preserve CIGAR interpretation, qualities, HPC lengths and 250-read cap |
| `alignment.rs` quality-error estimation and posterior analysis | `Consensus.Posterior` | reference-required | Match log-domain calculations, trimming, masking and low-quality criteria |
| `alignment.rs::plan_pileup_batches` | `Consensus.BatchPlanner` | support | Memory scheduling may vary only if results and deterministic order remain identical |
| `alignment.rs::merge_similar_consensuses` | `Refinement.Merge` | reference-required | Must be characterized branch by branch before coding |
| `chimera.rs` | `Refinement.Chimeras` | reference-required | Minimap2/CIGAR semantics, depth ratio, parent similarity and coverage thresholds |
| `alignment.rs` ASV mapping and EM | `Quantification.Mapping`, `Quantification.EM` | reference-required | Map candidate construction, likelihood/weights, convergence, depth rounding |
| `alignment.rs` pooled sample quantification | `Quantification.Pooled` | reference-required | Shared discovery, independent per-file abundance; preserve sort-index invariant |
| `classify.rs` | `Taxonomy.AlignmentClassifier`, `Taxonomy.EM` | reference-required | Separate EM from ASV EM; preserve best-hit and threshold behavior |
| `sintax.rs` | `Taxonomy.Sintax` | reference-required | Custom xorshift, sampling with replacement, 12-mers, denominator and tie behavior |
| `taxonomy.rs` | `Taxonomy.Database`, `Taxonomy.Output` | reference-required | EMU/SILVA/GG2/UNITE parsers require fixture coverage |
| `databases.rs`, `download.rs` | `Taxonomy.Registry`, `Taxonomy.Downloads` | support | Prefer Julia `Downloads`/archive packages; checksums and atomic installation should be added |
| `merge.rs` | `Export` | reference-required | Exact RC-canonical hashing and optional terminal fuzzy merge |
| `utils.rs` | functions in owning modules; `Diagnostics.Memory` | support | Avoid a generic utility dumping ground |
| `mapping.rs`, `map_processing.rs` | none initially | inactive-upstream | Commented out in `lib.rs`; retain only as historical reference |
| `polishing_mod.rs`, `polishing/` | none initially | inactive-upstream | Current POA is in `alignment.rs` using `spoars` |
| `graph.rs`, `unitig.rs` | none initially | inactive-upstream | Large historical assembly subsystem, not part of 0.8.0 binary |

## External dependency mapping

| Rust dependency | Upstream role | Julia reference candidate | Decision for reference phase |
|---|---|---|---|
| `needletail` + `flate2` | FASTA/FASTQ and compressed input | `FASTX.jl`, `CodecZlib.jl`; `CodecBzip2.jl` if bz2 parity is required | Use existing packages; test multiline FASTA, FASTQ qualities, gzip and header fidelity |
| `bio-seq` | packed DNA and reverse complements | `BioSequences.jl` plus explicit compact marker structs | Use existing package for DNA; benchmark `LongDNA{2}`/`LongDNA{4}` versus custom words |
| `fxhash` | fast non-cryptographic maps/sets | Base `Dict`; evaluate `DataStructures.RobinDict` and `SwissDict` | Select by measured lookup/iteration cost; never depend on hash iteration order |
| `smallvec` | inline short vectors | `StaticArrays.jl`, `SmallCollections.jl`, or fixed tuples | Choose per concrete maximum size after type/allocation tests |
| `fastbloom` | candidate prefilter | `BloomFilters.jl` only if bit/hash semantics can match; otherwise small internal implementation | Exact false-positive/hash behavior can affect retained candidates, so validate carefully |
| `crossbeam-channel`, `spmc` | bounded streaming pipelines | `Channel`, tasks, thread-local buffers | Julia-native; benchmark backpressure and GC behavior |
| `rayon` | data parallelism | `Threads.@threads`, `Threads.@spawn`, perhaps `ChunkSplitters.jl` | Start with Base threads and deterministic per-thread accumulators |
| `fishers_exact` | two-sided Fisher test | `HypothesisTests.jl` or an exact local kernel | Validate tail definition exhaustively; package API convenience alone is insufficient |
| `statrs` | binomial CDF | `Distributions.jl` | Already used by local DADA2.jl; compare boundary/rounding cases |
| `spoars` | banded/exact overlap POA | no proven drop-in Julia package | Fidelity spike first; native compatibility backend may be required temporarily |
| `minimap2` Rust binding | mapping, CIGAR, identity, database index | isolated `Minimap2Backend`; native binding/artifact initially | Pure Julia replacement only after differential equivalence is demonstrated |
| `clap`, `flexi_logger`, `log` | CLI and logging | `ArgParse.jl` or `Comonicon.jl`; `Logging`, `LoggingExtras.jl` | CLI choice does not enter core algorithms; structured stage metadata required |
| `serde`, `bincode` | serialization, partly historical | `Serialization` is not a stable interchange format; use versioned JLD2/Arrow only if checkpoints are implemented | Deferred until checkpoint schema is designed |
| `memory-stats`, `sysinfo` | diagnostics | Base/OS probes, `MemoryInspector.jl` for development | Keep optional; Linux peak RSS also captured externally |
| `clap-markdown` | generated CLI documentation | Documenter examples generated from Julia CLI schema | Support only |

Official package documentation confirms that FASTX records are byte-oriented and should be converted to a typed BioSequence when validation is needed; this supports the proposed parser/packed-sequence boundary. DataStructures supplies Robin Hood and Swiss-table dictionaries, but they remain benchmark candidates rather than assumed improvements. KernelAbstractions can target multiple accelerators, but its own documentation notes GPU-style semantics can impose extra CPU work; it is therefore deferred from the reference mode.

## Local package reuse assessment

| Local package | Relevant capability | Reference-mode decision |
|---|---|---|
| `DADA2.jl` | BioSequences/CodecZlib/Distributions dependency stack; denoising test patterns | Reuse dependency knowledge and fixtures, not DADA2's different inference algorithm |
| `ProfileAlignments.jl` | native Julia pair/profile alignment, k-mer distances, progressive alignment | Useful oracle/experimental backend; not equivalent to SPOA overlap POA or minimap2 |
| `AASV.jl` | classification and sequence utilities | Reuse only individually validated generic parsers/utilities; avoid circular package dependency |
| Elbrecht outputs | downstream ASV analysis tables | Use as an integration fixture after core parity; it is data, not a package dependency |

## Fidelity risk register

| Risk | Why it can change biology | Required evidence before acceptance |
|---|---|---|
| POA replacement | changes indels, endpoints and single-base consensus | exhaustive small POA fixtures plus real cluster differential tests |
| minimap2 replacement | changes candidate mappings, CIGARs, chimera parents and EM edges | mapping-level parity on synthetic substitutions/indels/RC/end clipping |
| statistical-tail mismatch | changes which allele pairs become SNPmers | exhaustive small contingency/binomial tables and boundary tests |
| hash/parallel order | changes greedy representatives and tied assignments | byte-identical outputs across runs and thread counts |
| quality representation | changes filtering and posterior consensus | golden quality-bin, identity and posterior cases |
| sparse EM refactor | changes summation order/convergence | trace-level comparison of every iteration under specified tolerances |
| streaming/chunking | can alter global counts and greedy order | invariance tests over chunk sizes and input-file partitions |
