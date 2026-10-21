# Upstream oracle

`rust_upstream` compiles the pinned savont 0.8.0 modules (`cli`, `constants`,
`utils`, `types`, `seeding`, `seq_parse`, `kmer_comp`, `asv_cluster`, `alignment`)
unchanged, with minimap2 2.30 (crate 0.1.31) and spoars 0.1.3 as pinned upstream, through
`#[path]` includes of `reference/savont/src`, and prints their outputs as
tab-separated records. Only `databases` is replaced by a stub holding the download
keywords that `cli.rs` needs.

## Build

1. Clone the pinned upstream checkout (ignored by git):

   ```sh
   git clone https://github.com/bluenote-1577/savont reference/savont
   git -C reference/savont checkout 5a5a1ce9053a867abc09e991f61b8ec7fb4fc6ef
   ```

2. Fetch the crates and patch the minimap2 wrapper (its index-path handling uses a
   Unix-only API; alignment code is untouched):

   ```sh
   (cd oracle/rust_upstream && cargo fetch)
   sh oracle/vendor_minimap2.sh
   ```

   On Windows, `oracle/rust_upstream/.cargo/config.toml` links the zlib of
   MinGW (Strawberry Perl).

3. Build from `oracle/rust_upstream`:

   ```sh
   cargo build --release
   ```

   Without the MSVC linker on Windows, use the GNU toolchain:
   `cargo +stable-x86_64-pc-windows-gnu build --release`. Add `--offline` when
   the crates are already cached. `Cargo.lock` resolves the upstream dependency
   ranges to statrs 0.16.1, fishers_exact 1.0.1, minimap2 0.1.31 and spoars 0.1.3.

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
  sorted before printing because upstream collects them in thread order;
- `cluster <savont asv arguments>`: retained reads, primary clusters, greedy
  SNPmer clusters and final stage-3 clusters;
- `poa <fastq> <band|noband>`: spoars alignment of every record to the growing
  graph with savont's settings, and the consensus. `SPOARS_FORCE_ISA=sse41`
  selects 128-bit vectors.
