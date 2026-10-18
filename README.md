# Savont.jl

Native Julia reimplementation of [savont](https://github.com/bluenote-1577/savont), pinned to upstream 0.8.0 at commit `5a5a1ce9053a867abc09e991f61b8ec7fb4fc6ef`.

The implementation advances through differential-validation gates. The current gate provides configuration, FASTA/FASTQ and gzip input, Cutadapt `rc` handling, sequence/HPC utilities, 2-bit k-mer encoding, the reversible minimap hash, FracMinHash seeds, minimizers, SNPmers, and blockmers. The end-to-end clustering pipeline is not yet released.

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

`test/rust_oracle.jl` optionally compares primitive outputs with the pinned Rust extraction when `SAVONT_RUST_ORACLE` points to the compiled oracle executable. The extraction avoids linking the unrelated minimap2 dependency, which fails in the upstream Windows build.

## Attribution

The original savont software is by Jim Shaw and contributors. The extracted reference routines retain the upstream MIT notice in `LICENSE` and identify their source commit.

Shaw J, Riisgaard-Jensen M, Andersen KS, Kirkegaard RH, Dueholm MKD, Li H. *Sensitive long-read amplicon sequence variant recovery with savont*. bioRxiv (2026). DOI: 10.64898/2026.05.26.727271.
