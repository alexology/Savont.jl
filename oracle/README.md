# Upstream oracle

`rust_upstream` compiles the pinned savont 0.8.0 modules (`chimera`, `cli`, `constants`,
`utils`, `types`, `seeding`, `seq_parse`, `kmer_comp`, `asv_cluster`, `alignment`,
`taxonomy`, `databases`, `classify`, `sintax`, `merge`) unchanged, with minimap2 2.30
(crate 0.1.31) and spoars 0.1.3 as pinned upstream, through `#[path]` includes of
`reference/savont/src`, and prints their outputs as tab-separated records. Only
`main.rs` (which cannot be included) and `download.rs` are left out.

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
  selects 128-bit vectors;
- `consensus <savont asv arguments>`: stages 1–3 as `cluster`, then
  `alignment::align_and_consensus` (minimap2 orientation and spoars POA): one line
  per initial consensus with cluster index, depth, read indices and sequence;
- `seedmap <fastx>`: maps every record after the first to the first with
  minimap2 configured as in `align_and_consensus` (strand, coordinates, mapping
  quality, CIGAR with soft clips), to check Minimap2.jl in that configuration;
- `polish <savont asv arguments>`: stages 1–4 as upstream `main.rs` (pileups of
  the deepest 10% for the quality error rates, then pileup analysis batch by
  batch and the low-quality partition). Prints the error rates, a per-position
  pileup summary (`depth:matches:deletions:insertions:quality sum:inserted
  quality sum:median HP length`) and each polished consensus (low-quality flag,
  HP lengths, masked sequence, low-quality positions, decompressed sequence);
- `merge <savont asv arguments>`: stages 1–4, decompression, then
  `alignment::merge_similar_consensuses` (stage 5); one line per merged
  consensus with index, depth, read indices and decompressed sequence;
- `chimera <savont asv arguments>`: stages 1–5, then `chimera::detect_chimeras`
  and `filter_chimeras` (stage 6): chimera score of every consensus, the
  detected combinations (sorted) and the remaining consensuses;
- `em <savont asv arguments>`: stages 1–6, then
  `alignment::refine_asv_depths_with_em` (stage 7) and the final sort: the ASVs
  with depth and assignment counts, and the sorted rows of
  `read_to_asv_mappings.tsv`;
- `asv <savont asv arguments>`: the whole pipeline as `main.rs::run_cluster`
  (stage 7b with `--pooled-samples`), printing `final_asvs.fasta`,
  `feature-table.tsv` and `final_assignments.tsv`, or `stopped` where `main.rs`
  returns early (`--skip-chimera-detection`). `write_feature_table` is private to
  `main.rs` and is reproduced in the oracle.
- `taxonomy <savont classify|sintax|export arguments>`: the subcommand as
  `main.rs` runs it (thread pool and database loading, without the log file),
  then every file of its output directory as `<file name>	<line>`, in name order.
  `test/oracle_cases.jl` copies the database to a temporary directory first,
  because classify caches its minimap2 index (`<fasta>.mmi`) next to the FASTA, and
  passes `-t 1` so that sintax ties between references follow database order.
