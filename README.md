# Savont.jl

Native Julia reimplementation of [savont](https://github.com/bluenote-1577/savont), pinned to upstream 0.8.0 at commit `5a5a1ce9053a867abc09e991f61b8ec7fb4fc6ef`.

The implementation advances through differential-validation gates. Implemented and checked against the compiled upstream code:

- **Phase 3A, primitives:** configuration, streaming FASTA/FASTQ and gzip input, Cutadapt `rc` handling, sequence/HPC utilities, 2-bit k-mer encoding, the reversible minimap hash, FracMinHash seeds, minimizers, split k-mers and blockmers;
- **Phase 3B, marker discovery:** strand-separated split k-mer and blockmer counting (`count_split_kmers`), the high-frequency k-mer threshold and biallelic SNPmer/blockmer inference with upstream's binomial and Fisher tests (`find_snpmers`, `find_blockmers`, `discover_markers`). The tests are exact ports of statrs 0.16.1 and fishers_exact 1.0.1 (`binomial_test`, `fishers_exact`). This reference path is single-threaded; Bloom prefiltering (`bloom_filter_size > 0`, off by default upstream) is not implemented.

Read construction (`TwinRead`), clustering, consensus, chimera detection, quantification and taxonomy are not implemented yet, so there is no end-to-end pipeline.

Known upstream bugs that Savont.jl corrects, with an upstream-compatible mode kept for testing, are listed in [`docs/deviations.md`](docs/deviations.md).

The Phase-1 design is recorded in:

- [`docs/architecture.md`](docs/architecture.md)
- [`docs/algorithm_mapping.md`](docs/algorithm_mapping.md)
- [`docs/implementation_plan.md`](docs/implementation_plan.md)
- [`docs/validation_plan.md`](docs/validation_plan.md)

Run the Julia tests with:

```julia
using Pkg
Pkg.test()
```

`test/rust_oracle.jl` compares the Julia results with outputs of the upstream modules compiled unchanged (`oracle/rust_upstream`, see [`oracle/README.md`](oracle/README.md)) on deterministic synthetic inputs in `test/data`: about 15,000 checks, including random sequences with ambiguous bases, three k-mer settings and three marker-discovery configurations. Without the oracle, the stored outputs in `test/reference` are used; with `SAVONT_RUST_ORACLE` pointing to the compiled oracle, it is run live.

## Attribution

The original savont software is by Jim Shaw and contributors. The extracted reference routines retain the upstream MIT notice in `LICENSE` and identify their source commit.

Shaw J, Riisgaard-Jensen M, Andersen KS, Kirkegaard RH, Dueholm MKD, Li H. *Sensitive long-read amplicon sequence variant recovery with savont*. bioRxiv (2026). DOI: 10.64898/2026.05.26.727271.
