# Savont.jl

Native Julia reimplementation of [savont](https://github.com/bluenote-1577/savont), pinned to upstream 0.8.0 at commit `5a5a1ce9053a867abc09e991f61b8ec7fb4fc6ef`: amplicon sequence variants (ASVs) from long reads (ONT, PacBio HiFi). All five subcommands are implemented: the `savont asv` pipeline, taxonomic classification (`classify`, `sintax`), `export` and `download`; on the test data they give the same output files as the compiled upstream code.

This package is part of an experimental effort to use AI to port
bioinformatics and quantitative ecology packages to Julia, with the goal of
creating a unified environment for high performance scientific computing.

## Install from GitHub

This package is not registered in the Julia General registry. Install it
from GitHub together with its unregistered dependency:

```julia
using Pkg
Pkg.add(url="https://github.com/alexology/Minimap2.jl")
Pkg.add(url="https://github.com/alexology/Savont.jl")
```

## Usage

```julia
using Savont
config = SavontConfig(min_read_length=1100, max_read_length=2000)   # upstream defaults
result = savont_asv(["sample.fastq.gz"], config; output_dir="savont-out")
result.asvs        # final ASVs (ConsensusSequence), by decreasing depth
```

Taxonomy and export of one or more runs:

```julia
download_database("databases", "emu-1")                       # also silva-138.2, greengenes2-2024.09, unite-10.0
savont_classify("savont-out", "databases/emu-1")              # species_abundance.tsv, genus_abundance.tsv, asv_mappings.tsv
savont_sintax("savont-out", "databases/emu-1"; output_dir="savont-out/sintax")   # genus level, k-mer bootstrap
savont_export(["run1-out", "run2-out"], "merged")             # QIIME2 tables, ASVs identified by sequence
```

From the command line, with upstream's subcommands and options (`--help` lists them; the number of threads is set by `julia -t`):

```sh
julia -t 2 --project=Savont.jl Savont.jl/bin/savont.jl asv sample1.fastq.gz sample2.fastq.gz --pooled-samples -m 2700 -M 4000 -o savont-out
julia -t 2 --project=Savont.jl Savont.jl/bin/savont.jl classify -i savont-out -d databases/emu-1
julia --project=Savont.jl Savont.jl/bin/savont.jl export -i run1-out run2-out -o merged
```

As upstream, reads must be primer/adapter trimmed, and `min_read_length`/`max_read_length` must bracket the amplicon length. Presets: `apply_preset(config, :hifi)` and `apply_preset(config, :rrna_operon)`. With several input files, `pooled_samples=true` adds per-sample depths. `output_dir` receives `final_asvs.fasta`, `feature-table.tsv` (QIIME2 format) and `final_assignments.tsv`, in upstream's formats.

Start Julia with threads (`julia -t auto`): read processing, marker clustering, consensus, polishing, merging, chimera detection and quantification run in parallel, and results do not depend on the number of threads. Primary clustering is sequential, as upstream's greedy algorithm. Memory: a retained read takes about 6 KB (2-bit bases and 4-bit quality bins, as upstream); a whole ONT sample of 365,000 reads of 2.8 kb (2700–4000 bp window) ran with 2 threads at a peak of 4.4 GB resident memory (2.3 GB of live data) in 30 minutes.

## Stages

1. **Markers** (`discover_markers`): strand-separated split k-mer and blockmer counts and biallelic SNPmer/blockmer inference with upstream's binomial and Fisher tests (exact ports of statrs 0.16.1 and fishers_exact 1.0.1). Bloom prefiltering (`bloom_filter_size > 0`, off by default upstream) is not implemented.
2. **Retained reads and primary clustering** (`build_twin_reads`, `cluster_reads_by_kmers`): length, identity and repetitiveness filters, open-syncmer minimizers, SNPmers, blockmers, MinHash LSH signatures; greedy k-mer clustering.
3. **Polymorphic-marker clustering** (`cluster_reads_by_snpmers`) with iterative consensus reclustering.
4. **Consensus and polishing**: per cluster, reads oriented with [Minimap2.jl](https://github.com/alexology/Minimap2.jl) (native port of minimap2 2.30, configured as upstream's minimap2 crate) and combined by partial order alignment (`initial_consensus`; a port of spoars 0.1.3 with its banded SIMD kernels emulated lane by lane); pileups, per-quality error rates, posterior analysis, masking and low-quality consensuses (`polish_consensuses!`).
5. **Merging** of similar consensuses (`merge_similar_consensuses`).
6. **Chimera detection** (`detect_chimeras`, `filter_chimeras`).
7. **Quantification** (`refine_depths!`): reads assigned to ASVs by SNPmers and alignment (or by minimap2 in low-polymorphism mode), EM over equivalence classes; per-sample depths with pooled samples (`per_sample_depths`).

## Taxonomy and export

- **Databases** (`load_database`, `download_database`): EMU, SILVA 138.2, GreenGenes2 2024.09 and UNITE 10.0 in upstream's layouts, typed by the `.savont_db` marker or the directory name. Downloads use Julia's `Downloads` and `Tar` instead of `wget`, `tar` and `gzip`; the UNITE FASTA is rewritten with short names as upstream.
- **`savont_classify`**: ASVs aligned to the database with Minimap2.jl (`map-ont`, as upstream's minimap2 crate), all alignments with the NM of the primary one kept, EM over taxa for multi-mapped ASVs, lineages cut at the rank their identity reaches (species ≥ 99%, genus ≥ 94.5%, then 86.5/82/78.5/75%). Pooled runs give one column per sample. The index is built in memory at every call (upstream caches it as `<fasta>.mmi`), so very large databases such as SILVA need a lot of memory.
- **`savont_sintax`**: upstream's SINTAX (canonical 12-mers, 32 drawn per iteration by its xorshift generator, 100 iterations, bootstrap ≥ 0.8), genus level and above.
- **`savont_export`**: ASVs of several runs combined by a strand-independent sequence hash, with fuzzy merging of copies trimmed by up to 10 bases, lineages from `asv_mappings.tsv`, and `merged_feature_table.tsv`, `merged_rep_seqs.fasta`, `merged_asv_taxonomy.tsv`, `merged_taxon_counts.tsv`.

Classification and export run on several threads (alignments, SINTAX references); their results do not depend on the number of threads. On a synthetic database of EMU's size (50,000 references of 1,450 bases) with 200 ASVs and 4 threads, the output files equal upstream's; classify takes about 9 s after compilation (upstream 2.9 s, or 1.2 s with its cached index), mostly to build the Minimap2.jl index (7 s, peak memory 3.8 GB), and sintax 1.2 s (upstream 1.0 s).

## Deviations from upstream

Upstream bugs that Savont.jl corrects are listed, with their measured effects, in [`docs/deviations.md`](docs/deviations.md); `SavontConfig(upstream_quirks=true)` reproduces upstream exactly and is used for the comparisons. The ones that change results on real data are D6 (the pairwise similarities of chimera detection are never found upstream, so two of its rules are inactive; the correction, chosen as the default after review, activates rules that were never exercised upstream: on 100,000 real reads 28 ASVs instead of 25) and D9 (`final_assignments.tsv` reports every read as unassigned when the FASTQ headers carry tab-separated tags, as dorado's do). D11–D15 (2026-10-22) reject even sampling rates (they separate the strands), keep EM depths summing to the assigned reads (one read more on one ASV of 21 on real data), give reads with repeated names their own assignments, avoid stale tables after empty reruns and reject a zero total depth in classify. D1 also changes `export`: copies of an ASV trimmed at the 5' end are merged, as intended upstream. Where upstream output depends on hash-map order or thread timing (ties in abundance tables, SINTAX votes and scores, fuzzy-merge representatives), Savont.jl uses a fixed order, described at the end of `docs/deviations.md`.

## Validation

`test/rust_oracle.jl` compares the Julia results with outputs of the upstream modules compiled unchanged (`oracle/rust_upstream`, see [`oracle/README.md`](oracle/README.md)) on deterministic synthetic inputs in `test/data` (about 64,000 checks; 64,195 tests in all): primitives on random sequences with ambiguous bases, statistics, three marker-discovery configurations, and fourteen pipeline configurations (k-mer settings, blockmers, single strand, low polymorphism, reads with 1–2% substitutions, mixed haplotypes with homopolymer compression and masking, chimeras), compared stage by stage: every field of every retained read, primary, greedy and final clusters, partial order alignments (exact and banded, 256- and 128-bit vectors), initial consensuses, pileups position by position, error rates, polished and merged consensuses, chimera calls, ASV depths and read assignments, and the three output files byte for byte (also with pooled samples and with `--skip-chimera-detection`). For the taxonomy subcommands, synthetic databases in the four formats and ASV directories (`test/data/taxonomy`, generated by `tools/make_test_data.jl taxonomy`) reach every identity band, ties between identical references, EM redistribution, reverse strands, ambiguous and lower-case bases, keys missing from the taxonomy and pooled samples; the output files of 7 classify, 5 sintax and 4 export runs are identical to upstream's (rows of equal abundance compared as sets, see above). Downloads are tested offline on local archives. Without the oracle, the stored outputs in `test/reference` are used; with `SAVONT_RUST_ORACLE` pointing to the compiled oracle, it is run live.

The pipeline was also compared with the oracle on real ONT reads of a 2.8 kb amplicon: on 20,000 reads the three output files are identical (20 ASVs, 13,683 reads), and stages 1–4 are identical on 100,000 reads.

Minimap2.jl is not registered; `Project.toml` takes it from `../Minimap2.jl` (`[sources]`), so both repositories must be checked out side by side. Run the tests with:

```julia
using Pkg
Pkg.test()
```

The Phase-1 design is recorded in [`docs/architecture.md`](docs/architecture.md), [`docs/algorithm_mapping.md`](docs/algorithm_mapping.md), [`docs/implementation_plan.md`](docs/implementation_plan.md) and [`docs/validation_plan.md`](docs/validation_plan.md).

## Attribution

The original savont software is by Jim Shaw and contributors. The extracted reference routines retain the upstream MIT notice in `LICENSE` and identify their source commit.

Shaw J, Riisgaard-Jensen M, Andersen KS, Kirkegaard RH, Dueholm MKD, Li H. *Sensitive long-read amplicon sequence variant recovery with savont*. bioRxiv (2026). DOI: 10.64898/2026.05.26.727271.
