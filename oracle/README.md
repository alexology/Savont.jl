# Upstream oracle

`rust_upstream` compiles the pinned savont 0.8.0 modules (`cli`, `constants`,
`utils`, `types`, `seeding`, `seq_parse`, `kmer_comp`) unchanged, through
`#[path]` includes of `reference/savont/src`, and prints their outputs as
tab-separated records. Only `databases` is replaced by a stub holding the download
keywords that `cli.rs` needs. minimap2 and spoars are not linked, so the oracle
builds on Windows.

## Build

1. Clone the pinned upstream checkout (ignored by git):

   ```sh
   git clone https://github.com/bluenote-1577/savont reference/savont
   git -C reference/savont checkout 5a5a1ce9053a867abc09e991f61b8ec7fb4fc6ef
   ```

2. Build from `oracle/rust_upstream`:

   ```sh
   cargo build --release
   ```

   Without the MSVC linker on Windows, use the GNU toolchain:
   `cargo +stable-x86_64-pc-windows-gnu build --release`. Add `--offline` when
   the crates are already cached. `Cargo.lock` resolves the upstream dependency
   ranges to the cached versions (statrs 0.16.1, fishers_exact 1.0.1).

## Use

```sh
SAVONT_RUST_ORACLE=oracle/rust_upstream/target/release/savont-upstream-oracle \
    julia --project=. -e 'using Pkg; Pkg.test()'
```

compares against the live oracle instead of the stored outputs in
`test/reference`. Regenerate those with `tools/make_reference.jl` after changing
`test/oracle_cases.jl`, and the inputs in `test/data` with
`tools/make_test_data.jl`.

Commands:

- `primitives <fastx> <k> <c> <w> <l> <minimum_bq>`: per-record FracMinHash seeds,
  minimizers, split k-mers, blockmers, reverse complement, homopolymer
  compression, quality pooling and identity estimates;
- `stats`: binomial tail and Fisher exact p-values on a fixed grid;
- `markers <savont asv arguments>`: `read_to_split_kmers`,
  `get_blockmers_inplace_sort` and `get_snpmers_inplace_sort`. Blockmers are
  sorted before printing because upstream collects them in thread order.
