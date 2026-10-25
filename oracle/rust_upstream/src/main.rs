//! Differential oracle built from the pinned upstream savont 0.8.0 sources.
//! Source commit: 5a5a1ce9053a867abc09e991f61b8ec7fb4fc6ef.
//!
//! The upstream modules are compiled unchanged; this file only drives them and
//! prints tab-separated records for test/rust_oracle.jl. Driver code copied from
//! main.rs (which cannot be included) is marked as such.
#![allow(dead_code, unused_imports, unused_variables, unused_mut)]

#[path = "../../../reference/savont/src/alignment.rs"]
mod alignment;
#[path = "../../../reference/savont/src/asv_cluster.rs"]
mod asv_cluster;
#[path = "../../../reference/savont/src/chimera.rs"]
mod chimera;
#[path = "../../../reference/savont/src/cli.rs"]
mod cli;
#[path = "../../../reference/savont/src/constants.rs"]
mod constants;
#[path = "../../../reference/savont/src/kmer_comp.rs"]
mod kmer_comp;
#[path = "../../../reference/savont/src/seeding.rs"]
mod seeding;
#[path = "../../../reference/savont/src/seq_parse.rs"]
mod seq_parse;
#[path = "../../../reference/savont/src/types.rs"]
mod types;
#[path = "../../../reference/savont/src/utils.rs"]
mod utils;

#[path = "../../../reference/savont/src/classify.rs"]
mod classify;
#[path = "../../../reference/savont/src/databases.rs"]
mod databases;
#[path = "../../../reference/savont/src/merge.rs"]
mod merge;
#[path = "../../../reference/savont/src/sintax.rs"]
mod sintax;
#[path = "../../../reference/savont/src/taxonomy.rs"]
mod taxonomy;

use clap::Parser;

fn csv<T: ToString>(values: &[T]) -> String {
    values.iter().map(T::to_string).collect::<Vec<_>>().join(",")
}

/// `primitives <fastx> <k> <c> <w> <l> <minimum_bq>`: per-record primitive outputs.
fn primitives(args: &[String]) {
    let path = &args[0];
    let k: usize = args[1].parse().unwrap();
    let c: usize = args[2].parse().unwrap();
    let w: usize = args[3].parse().unwrap();
    let l: usize = args[4].parse().unwrap();
    let minimum_bq: u8 = args[5].parse().unwrap();
    let mut reader = needletail::parse_fastx_file(path).expect("valid path");
    let mut index = 0usize;
    while let Some(record) = reader.next() {
        index += 1;
        let record = record.expect("valid record");
        let seq = record.seq().to_vec();
        let qualities = record.qual().map(Vec::from);

        let mut hashes = Vec::new();
        seeding::fmh_seeds(&seq, &mut hashes, c, k);
        println!("fmh\t{index}\t{}", csv(&hashes));

        let mut kmers = Vec::new();
        let mut positions = Vec::new();
        seeding::fmh_seeds_positions(&seq, &mut kmers, &mut positions, c, k);
        println!("fmhp\t{index}\t{}\t{}", csv(&kmers), csv(&positions));

        let mut minimizers = Vec::new();
        let mut minimizer_positions = Vec::new();
        seeding::minimizer_seeds_positions(&seq, &mut minimizers, &mut minimizer_positions, w, k);
        println!("minimizer\t{index}\t{}\t{}", csv(&minimizers), csv(&minimizer_positions));

        let split = seeding::split_kmer_mid(seq.clone(), qualities.clone(), k, minimum_bq);
        println!("split\t{index}\t{}", csv(&split));

        let blocks = seeding::blockmer_kmers(seq.clone(), qualities.clone(), k, l, minimum_bq);
        let blocks = blocks
            .iter()
            .map(|x| format!("{}:{}", x.kmer, u8::from(x.is_forward)))
            .collect::<Vec<_>>();
        println!("blocks\t{index}\t{}", blocks.join(","));

        println!("revcomp\t{index}\t{}", String::from_utf8(utils::reverse_complement(&seq)).unwrap());
        let (hpc, lengths) = utils::homopolymer_compress(&seq, true);
        println!("hpc\t{index}\t{}\t{}", String::from_utf8(hpc).unwrap(), csv(&lengths));

        if let Some(q) = qualities.as_ref() {
            let (hpc, hpc_q, lengths) = utils::homopolymer_compress_with_quality(&seq, q, true);
            println!(
                "hpcq\t{index}\t{}\t{}\t{}",
                String::from_utf8(hpc).unwrap(),
                csv(&hpc_q),
                csv(&lengths)
            );
            println!("pool\t{index}\t{}", csv(&seq_parse::quality_pool(q.clone())));
            println!("identity\t{index}\t{:?}", seeding::estimate_sequence_identity(Some(q)).unwrap());
        }
    }
}

/// `stats`: binomial tail and Fisher exact p-values on a fixed grid.
fn stats() {
    for p in [0.025_f64, 0.05] {
        for n in (1..=60u64).chain([100, 250, 1000, 5000, 20000]) {
            let step = (n / 60).max(1);
            let mut succ = 0;
            while succ <= n {
                println!("binomial\t{n}\t{succ}\t{p:?}\t{:?}", utils::binomial_test(n, succ, p));
                succ += step;
            }
        }
    }
    let values = [0u32, 1, 2, 3, 5, 8, 13, 40, 120];
    for &a in &values {
        for &b in &values {
            for &c in &values {
                for &d in &values {
                    let table = [a, b, c, d];
                    let p = fishers_exact::fishers_exact(&table).unwrap().two_tail_pvalue;
                    println!("fisher\t{a}\t{b}\t{c}\t{d}\t{p:?}");
                }
            }
        }
    }
    for table in [[1000u32, 15, 980, 30], [3500, 400, 3400, 380], [52, 0, 61, 0], [10000, 9000, 1, 2]] {
        let p = fishers_exact::fishers_exact(&table).unwrap().two_tail_pvalue;
        println!("fisher\t{}\t{}\t{}\t{}\t{p:?}", table[0], table[1], table[2], table[3]);
    }
}

/// `markers <savont asv arguments>`: k-mer counting plus SNPmer and blockmer inference.
fn markers(args: &[String]) {
    let mut argv = vec!["savont-asv".to_string()];
    argv.extend_from_slice(args);
    let cluster_args = cli::ClusterArgs::parse_from(argv);
    let k = cluster_args.kmer_size;
    let l = cluster_args.blockmer_length;
    let (snpmer_counts, blockmer_counts) =
        seq_parse::read_to_split_kmers(k, l, cluster_args.threads, &cluster_args);

    let mut sorted = snpmer_counts.clone();
    sorted.sort();
    for (kmer, counts) in sorted.iter() {
        println!("count\t{kmer}\t{}\t{}", counts[0], counts[1]);
    }
    let mut sorted_blocks = blockmer_counts.clone();
    sorted_blocks.sort();
    for (kmer, counts) in sorted_blocks.iter() {
        println!("blockcount\t{kmer}\t{}\t{}", counts[0], counts[1]);
    }

    if cluster_args.use_blockmers {
        let mut info =
            kmer_comp::get_blockmers_inplace_sort(blockmer_counts, &snpmer_counts, k, l, &cluster_args);
        // Upstream pushes from worker threads; sort for a deterministic comparison.
        info.blockmer_info.sort();
        for b in info.blockmer_info.iter() {
            println!(
                "blockmer\t{}\t{}\t{}\t{}\t{}",
                b.anchor_kmer, b.blockmers[0].kmer, b.blockmers[1].kmer, b.counts[0], b.counts[1]
            );
        }
    }

    let info = kmer_comp::get_snpmers_inplace_sort(snpmer_counts, k, &cluster_args);
    println!("high_freq_threshold\t{}", info.high_freq_thresh);
    let mut high_freq = info.high_freq_kmers.iter().map(|x| x.to_u64()).collect::<Vec<_>>();
    high_freq.sort();
    for kmer in high_freq {
        println!("high_freq\t{kmer}");
    }
    for s in info.snpmer_info.iter() {
        println!(
            "snpmer\t{}\t{}\t{}\t{}\t{}",
            s.split_kmer, s.mid_bases[0], s.mid_bases[1], s.counts[0], s.counts[1]
        );
    }
}

/// `poa <fastq> <band|noband>`: spoars partial order alignment with savont's
/// settings (alignment.rs::generate_consensus_poa): every record is aligned to the
/// growing graph and added with its raw quality bytes as weights.
fn poa(args: &[String]) {
    use spoars::align::{AlignmentEngine, AlignmentType, BandConfig, Scoring, SimdEngine};
    use spoars::graph::Graph;
    let mut reader = needletail::parse_fastx_file(&args[0]).expect("valid path");
    let mut sequences = Vec::new();
    let mut qualities = Vec::new();
    while let Some(record) = reader.next() {
        let record = record.expect("valid record");
        sequences.push(record.seq().to_vec());
        qualities.push(record.qual().expect("FASTQ input").to_vec());
    }
    let scoring = Scoring::new(3, -8, -6, -6, 0, 0).unwrap();
    // As generate_consensus_poa: band base = maximum deviation from the mean length.
    let ref_len = sequences.iter().map(|s| s.len()).sum::<usize>() / sequences.len();
    let max_deviation = sequences
        .iter()
        .map(|s| (ref_len as isize - s.len() as isize).abs() as u32)
        .max()
        .unwrap_or(0);
    let mut engine = if args[1] == "noband" {
        SimdEngine::new(AlignmentType::Overlap, scoring)
    } else {
        SimdEngine::banded(AlignmentType::Overlap, scoring, BandConfig { base: max_deviation, frac: 0.1 })
    };
    let mut graph = Graph::new();
    for (index, sequence) in sequences.iter().enumerate() {
        let weights: Vec<u32> = qualities[index].iter().map(|&q| q as u32).collect();
        let (alignment, score) = engine.align(sequence, &graph);
        let pairs = alignment.iter().map(|(n, s)| format!("{n}:{s}")).collect::<Vec<_>>();
        println!("alignment\t{}\t{score}\t{}", index + 1, pairs.join(","));
        graph.add_alignment(&alignment, sequence, &weights).expect("add_alignment failed");
    }
    println!("consensus\t{}", graph.generate_consensus());
}

/// `seedmap <fastx>`: maps every record after the first to the first one with the
/// minimap2 crate configured as in alignment.rs::align_and_consensus.
fn seedmap(args: &[String]) {
    let mut reader = needletail::parse_fastx_file(&args[0]).expect("valid path");
    let mut sequences = Vec::new();
    while let Some(record) = reader.next() {
        sequences.push(record.expect("valid record").seq().to_vec());
    }
    let aligner = minimap2::Aligner::builder().map_ont().with_index_threads(1).with_cigar().with_seq(&sequences[0]).unwrap();
    for (i, seq) in sequences.iter().enumerate().skip(1) {
        match aligner.map(seq, true, false, None, None, None) {
            Err(e) => println!("map\t{i}\terror\t{e}"),
            Ok(mappings) => {
                for m in mappings.iter() {
                    let cigar = m.alignment.as_ref().and_then(|a| a.cigar_str.clone()).unwrap_or_default();
                    println!("map\t{i}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", if m.strand == minimap2::Strand::Reverse { '-' } else { '+' },
                             m.query_start, m.query_end, m.target_start, m.target_end, m.mapq, cigar);
                }
            }
        }
    }
}

fn optional_csv<T: ToString>(values: &[Option<T>]) -> String {
    values
        .iter()
        .map(|x| x.as_ref().map_or("-".to_string(), T::to_string))
        .collect::<Vec<_>>()
        .join(",")
}

/// `cluster <savont asv arguments>`: twin reads (stage 1.5), primary k-mer
/// clustering (stage 2) and polymorphic-marker clustering (stage 3).
fn cluster(args: &[String]) {
    stages(args, true);
}

/// `consensus <savont asv arguments>`: stages 1-3, then the initial POA consensus
/// of each cluster (alignment.rs::align_and_consensus, stage 4a).
fn consensus(args: &[String]) {
    let (cluster_args, twin_reads, final_clusters, _) = stages(args, false);
    let directory = std::env::temp_dir().join(format!("savont-oracle-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let consensuses = alignment::align_and_consensus(&twin_reads, final_clusters, &cluster_args, &directory);
    std::fs::remove_dir_all(&directory).unwrap();
    for c in consensuses.iter() {
        println!("consensus\t{}\t{}\t{}\t{}", c.id, c.depth, csv(&c.cluster), String::from_utf8(c.sequence.clone()).unwrap());
    }
}

/// `polish <savont asv arguments>`: stages 1-4 as main.rs (initial consensus,
/// pileups, quality error rates, posterior analysis and low-quality partition).
/// Prints the error rates, a per-position summary of every pileup and each
/// polished consensus before the low-quality consensuses are removed.
fn polish(args: &[String]) {
    polish_stage(args, true);
}

/// `merge <savont asv arguments>`: stages 1-4, decompression as main.rs, then
/// stage 5 (alignment::merge_similar_consensuses). Prints the merged consensuses.
fn merge(args: &[String]) {
    let (cluster_args, twin_reads, consensuses, low_quality, _) = polish_stage(args, false);
    let directory = std::env::temp_dir().join(format!("savont-oracle-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let merged = alignment::merge_similar_consensuses(&twin_reads, consensuses, low_quality, &cluster_args, &directory);
    std::fs::remove_dir_all(&directory).unwrap();
    for c in merged.iter() {
        println!(
            "merged	{}	{}	{}	{}",
            c.id,
            c.depth,
            csv(&c.cluster),
            String::from_utf8(c.decompressed_sequence.clone().unwrap()).unwrap()
        );
    }
}

/// `chimera <savont asv arguments>`: stages 1-5, then chimera detection and
/// filtering (stage 6). Prints every consensus with its chimera score, the
/// detected chimera combinations (sorted; upstream collects them in parallel) and
/// the consensuses that remain.
fn chimera_stage(args: &[String]) {
    let (cluster_args, twin_reads, consensuses, low_quality, _) = polish_stage(args, false);
    let directory = std::env::temp_dir().join(format!("savont-oracle-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let mut merged = alignment::merge_similar_consensuses(&twin_reads, consensuses, low_quality, &cluster_args, &directory);
    std::fs::remove_dir_all(&directory).unwrap();
    let chimeras = chimera::detect_chimeras(&mut merged, &cluster_args);
    for c in merged.iter() {
        println!("score	{}	{}", c.id, c.chimera_score.map_or("-".to_string(), |s| s.to_string()));
    }
    let mut rows: Vec<String> = chimeras
        .iter()
        .map(|c| {
            format!(
                "chimera	{}	{}	{}	{}	{}	{}	{:?}",
                merged[c.query_idx].id,
                merged[c.left_parent_idx].id,
                merged[c.right_parent_idx].id,
                c.left_match_len,
                c.right_match_len,
                c.query_len,
                c.coverage_fraction
            )
        })
        .collect();
    rows.sort();
    for row in rows {
        println!("{}", row);
    }
    for c in chimera::filter_chimeras(merged, &chimeras).iter() {
        println!("kept	{}", c.id);
    }
}

/// `em <savont asv arguments>`: stages 1-6 (stage 6 unless --skip-chimera-detection),
/// then stage 7 (alignment::refine_asv_depths_with_em) and main.rs's final sort.
/// Prints the final ASVs in order and the sorted rows of read_to_asv_mappings.tsv.
fn em(args: &[String]) {
    let (cluster_args, twin_reads, consensuses, low_quality, kmer_info) = polish_stage(args, false);
    let directory = std::env::temp_dir().join(format!("savont-oracle-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let mut merged = alignment::merge_similar_consensuses(&twin_reads, consensuses, low_quality, &cluster_args, &directory);
    if !cluster_args.skip_chimera_detection {
        let chimeras = chimera::detect_chimeras(&mut merged, &cluster_args);
        merged = chimera::filter_chimeras(merged, &chimeras);
    }
    alignment::refine_asv_depths_with_em(&twin_reads, &mut merged, &kmer_info, &cluster_args, &directory);
    merged.sort_by(|a, b| b.depth.cmp(&a.depth).then_with(|| a.sequence.cmp(&b.sequence)).then_with(|| a.id.cmp(&b.id)));
    let mappings = std::fs::read_to_string(directory.join("read_to_asv_mappings.tsv")).unwrap_or_default();
    std::fs::remove_dir_all(&directory).unwrap();
    for c in merged.iter() {
        println!(
            "asv	{}	{}	{}	{}	{}	{}",
            c.id,
            c.depth,
            c.unambig_best_read_map_count.map_or("-".to_string(), |x| x.to_string()),
            c.ambig_read_map_count.map_or("-".to_string(), |x| x.to_string()),
            c.num_map_leq_10nm.map_or("-".to_string(), |x| x.to_string()),
            String::from_utf8(c.decompressed_sequence.clone().unwrap()).unwrap()
        );
    }
    let mut rows: Vec<&str> = mappings.lines().collect();
    rows.sort();
    for row in rows {
        println!("mapping	{}", row);
    }
}

/// `asv <savont asv arguments>`: the whole pipeline as main.rs::run_cluster from
/// stage 5 on (stage 7b with --pooled-samples), writing final_asvs.fasta,
/// feature-table.tsv and final_assignments.tsv to a temporary directory and
/// printing them (`fasta`, `table`, `assignment` lines). Prints `stopped` where
/// main.rs returns early (--skip-chimera-detection).
fn asv(args: &[String]) {
    let (cluster_args, twin_reads, consensuses, low_quality, kmer_info) = polish_stage(args, false);
    let directory = std::env::temp_dir().join(format!("savont-oracle-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let mut consensuses = alignment::merge_similar_consensuses(&twin_reads, consensuses, low_quality, &cluster_args, &directory);
    if cluster_args.skip_chimera_detection {
        std::fs::remove_dir_all(&directory).unwrap();
        println!("stopped");
        return;
    }
    let chimeras = chimera::detect_chimeras(&mut consensuses, &cluster_args);
    let mut consensuses = chimera::filter_chimeras(consensuses, &chimeras);
    alignment::refine_asv_depths_with_em(&twin_reads, &mut consensuses, &kmer_info, &cluster_args, &directory);
    consensuses.sort_by(|a, b| b.depth.cmp(&a.depth).then_with(|| a.sequence.cmp(&b.sequence)).then_with(|| a.id.cmp(&b.id)));
    let debug_to_final: std::collections::HashMap<usize, usize> =
        consensuses.iter().enumerate().map(|(i, c)| (c.id, i)).collect();
    let em_fasta = directory.join("final_asvs_for_em.fasta");
    alignment::write_consensus_fasta(&consensuses, &em_fasta, "em_refinement").unwrap();
    let sample_names: Vec<String> = cluster_args
        .input_files
        .iter()
        .map(|f| std::path::Path::new(f).file_stem().and_then(|s| s.to_str()).unwrap_or("sample").to_string())
        .collect();
    if cluster_args.pooled_samples && cluster_args.input_files.len() > 1 {
        let per_sample = alignment::compute_per_sample_depths(
            &twin_reads,
            cluster_args.input_files.len(),
            &consensuses,
            &kmer_info,
            &cluster_args,
            &em_fasta,
        );
        for (i, c) in consensuses.iter_mut().enumerate() {
            c.per_sample_depths = per_sample[i].clone();
        }
    }
    let fasta = directory.join("final_asvs.fasta");
    alignment::write_consensus_fasta(&consensuses, &fasta, "final").unwrap();
    // As main.rs::write_feature_table.
    let names: Vec<&str> = if cluster_args.pooled_samples && cluster_args.input_files.len() > 1 {
        sample_names.iter().map(|s| s.as_str()).collect()
    } else {
        vec![sample_names[0].as_str()]
    };
    let mut table = format!("#OTU ID	{}
", names.join("	"));
    for (i, c) in consensuses.iter().enumerate() {
        if c.per_sample_depths.is_empty() {
            let depth = c.depth + c.appended_depth;
            table.push_str(&format!("final_asv_{}_depth_{}	{}
", i, depth, depth));
        } else {
            let depths: Vec<String> = c.per_sample_depths.iter().map(|d| d.to_string()).collect();
            table.push_str(&format!("final_asv_{}_depth_{}	{}
", i, depths.join("-"), depths.join("	")));
        }
    }
    consensuses.iter_mut().enumerate().for_each(|(i, c)| c.id = i);
    let assignments = directory.join("final_assignments.tsv");
    alignment::write_final_assignments(
        &twin_reads,
        &consensuses,
        &directory.join("read_to_asv_mappings.tsv"),
        &debug_to_final,
        &assignments,
    )
    .unwrap();
    let fasta_text = std::fs::read_to_string(&fasta).unwrap();
    let assignment_text = std::fs::read_to_string(&assignments).unwrap();
    std::fs::remove_dir_all(&directory).unwrap();
    for line in fasta_text.lines() {
        println!("fasta	{}", line);
    }
    for line in table.lines() {
        println!("table	{}", line);
    }
    for line in assignment_text.lines() {
        println!("assignment	{}", line);
    }
}

/// Stages 1-4 up to the decompressed high- and low-quality consensuses.
fn polish_stage(
    args: &[String],
    verbose: bool,
) -> (
    cli::ClusterArgs,
    Vec<types::TwinRead>,
    Vec<types::ConsensusSequence>,
    Vec<types::ConsensusSequence>,
    types::KmerGlobalInfo,
) {
    let (cluster_args, twin_reads, final_clusters, kmer_info) = stages(args, false);
    let directory = std::env::temp_dir().join(format!("savont-oracle-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let mut consensuses = alignment::align_and_consensus(&twin_reads, final_clusters, &cluster_args, &directory);
    let batches = alignment::plan_pileup_batches(&consensuses, &cluster_args);
    let top_indices = alignment::top_cluster_indices(&consensuses, 0.1);
    let quality_error_map = {
        let top = alignment::generate_consensus_pileups(&twin_reads, &mut consensuses, &cluster_args, &top_indices);
        alignment::estimate_quality_error_rates(&top)
    };
    let mut qualities: Vec<u8> = quality_error_map.keys().cloned().collect();
    qualities.sort();
    for q in qualities.iter().filter(|_| verbose) {
        println!("rate	{}	{:?}", q, quality_error_map[q]);
    }
    for batch in &batches {
        let pileups = alignment::generate_consensus_pileups(&twin_reads, &mut consensuses, &cluster_args, batch);
        for (local, &index) in batch.iter().enumerate().filter(|_| verbose) {
            let summary: Vec<String> = pileups[local]
                .iter()
                .map(|p| {
                    let (mut matches, mut deletions, mut insertions, mut quality, mut inserted_quality) = (0, 0, 0, 0u64, 0u64);
                    for b in &p.bases {
                        match b {
                            alignment::PileupBase::Base(base, q, _) => {
                                matches += usize::from(*base == p.ref_base);
                                quality += *q as u64;
                            }
                            alignment::PileupBase::Deletion => deletions += 1,
                            alignment::PileupBase::Insertion(data) => {
                                insertions += 1;
                                inserted_quality += data.iter().map(|x| x.1 as u64).sum::<u64>();
                            }
                        }
                    }
                    format!("{}:{}:{}:{}:{}:{}:{}", p.depth(), matches, deletions, insertions, quality, inserted_quality, p.ref_hp_length)
                })
                .collect();
            println!("pileup	{}	{}", consensuses[index].id, summary.join(","));
        }
        alignment::analyze_pileup_consensuses(pileups, batch, &mut consensuses, &quality_error_map, &cluster_args);
    }
    let mut low_quality = alignment::collect_low_quality_consensuses(&mut consensuses, &twin_reads, &cluster_args, &directory);
    std::fs::remove_dir_all(&directory).unwrap();
    // As main.rs: decompress before merging.
    for c in consensuses.iter_mut().chain(low_quality.iter_mut()) {
        c.decompress();
    }
    if verbose {
        let low_ids: Vec<usize> = low_quality.iter().map(|c| c.id).collect();
        let mut all: Vec<&types::ConsensusSequence> = consensuses.iter().chain(low_quality.iter()).collect();
        all.sort_by(|a, b| b.depth.cmp(&a.depth).then_with(|| a.id.cmp(&b.id)));
        for c in all {
            println!(
                "polished	{}	{}	{}	{}	{}	{}	{}",
                c.id,
                c.depth,
                u8::from(low_ids.contains(&c.id)),
                csv(&c.hp_lengths),
                String::from_utf8(c.sequence.clone()).unwrap(),
                csv(&c.low_quality_positions),
                String::from_utf8(c.decompressed_sequence.clone().unwrap()).unwrap()
            );
        }
    }
    (cluster_args, twin_reads, consensuses, low_quality, kmer_info)
}

/// Stages 1-3; with `verbose`, prints retained reads and clusters.
fn stages(
    args: &[String],
    verbose: bool,
) -> (cli::ClusterArgs, Vec<types::TwinRead>, Vec<Vec<usize>>, types::KmerGlobalInfo) {
    let mut argv = vec!["savont-asv".to_string()];
    argv.extend_from_slice(args);
    let mut cluster_args = cli::ClusterArgs::parse_from(argv);
    let k = cluster_args.kmer_size;
    let l = cluster_args.blockmer_length;

    // As main.rs::get_kmers_and_snpmers.
    let (big_snpmer_map, big_blockmer_map) =
        seq_parse::read_to_split_kmers(k, l, cluster_args.threads, &cluster_args);
    let mut blockmer_info = if cluster_args.use_blockmers {
        kmer_comp::get_blockmers_inplace_sort(big_blockmer_map, &big_snpmer_map, k, l, &cluster_args)
    } else {
        types::BlockmerGlobalInfo::default()
    };
    let mut kmer_info = kmer_comp::get_snpmers_inplace_sort(big_snpmer_map, k, &cluster_args);

    // As main.rs::get_twin_reads_from_kmer_info.
    let mut twin_reads =
        kmer_comp::twin_reads_from_snpmers(&mut kmer_info, &mut blockmer_info, &cluster_args);
    twin_reads.sort_by(|a, b| {
        b.est_id
            .unwrap_or(100.0)
            .partial_cmp(&a.est_id.unwrap_or(100.0))
            .unwrap()
    });
    let without = twin_reads.iter().filter(|tr| tr.snpmers_vec().is_empty()).count();
    let fraction_without = without as f64 / twin_reads.len() as f64;
    if verbose {
        println!("low_polymorphism\t{}", u8::from(fraction_without > 0.75));
    }
    // As main.rs: auto-enable low-polymorphism mode.
    if fraction_without > 0.75 && !cluster_args.low_polymorphism {
        cluster_args.low_polymorphism = true;
    }

    for (index, read) in twin_reads.iter().enumerate().filter(|_| verbose) {
        let qualities = read
            .qual_seq
            .as_ref()
            .map_or("-".to_string(), |q| csv(&q.iter().map(|x| x as u8).collect::<Vec<_>>()));
        let est_id = read.est_id.map_or("-".to_string(), |x| format!("{x:?}"));
        let minimizer_kmers = read.minimizer_kmers.iter().map(|x| x.to_u64()).collect::<Vec<_>>();
        let snpmer_kmers = read.snpmer_kmers.iter().map(|x| x.to_u64()).collect::<Vec<_>>();
        let canonical = read.blockmer_canonical.iter().map(|&x| u8::from(x)).collect::<Vec<_>>();
        println!(
            "read\t{index}\t{}\t{}\t{}\t{}\t{est_id}\t{}\t{qualities}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            read.id,
            read.base_id,
            read.file_idx,
            read.base_length,
            read.dna_seq,
            csv(&read.minimizer_positions),
            csv(&minimizer_kmers),
            csv(&read.snpmer_positions),
            csv(&snpmer_kmers),
            csv(&read.blockmer_positions),
            csv(&canonical),
            optional_csv(&read.lsh_signatures)
        );
    }

    let directory = std::env::temp_dir().join(format!("savont-oracle-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let clusters = asv_cluster::cluster_reads_by_kmers(&twin_reads, &cluster_args, &directory);
    for members in clusters.iter().filter(|_| verbose) {
        println!("cluster\t{}", csv(members));
    }
    let final_clusters =
        asv_cluster::cluster_reads_by_snpmers(&twin_reads, &clusters, &cluster_args, &directory);
    // Greedy SNPmer clusters before reclustering, as written by upstream.
    let greedy = directory.join("snpmer_clusters_before_reclust2.5.tsv");
    if let (true, Ok(text)) = (verbose, std::fs::read_to_string(&greedy)) {
        for line in text.lines().skip(1) {
            let fields = line.split('\t').collect::<Vec<_>>();
            println!("greedy\t{}\t{}", fields[0], fields[4]);
        }
    }
    std::fs::remove_dir_all(&directory).unwrap();
    for members in final_clusters.iter().filter(|_| verbose) {
        println!("final\t{}", csv(members));
    }
    (cluster_args, twin_reads, final_clusters, kmer_info)
}

/// `taxonomy <savont classify|sintax|export arguments>`: runs the subcommand as
/// `main.rs` does (thread pool, database loading; logging omitted) and prints
/// every file of its output directory as `<file name>	<line>`, files in name order.
fn taxonomy_command(args: &[String]) {
    let cli = cli::Cli::parse_from(std::iter::once("savont".to_string()).chain(args.iter().cloned()));
    let pool = |threads: usize| {
        rayon::ThreadPoolBuilder::new().num_threads(threads).stack_size(16 * 1024 * 1024).build_global().unwrap()
    };
    let output_dir = match &cli.command {
        cli::Commands::Classify(a) => {
            pool(a.threads);
            let output = a.output_dir.clone().unwrap_or(a.input_dir.clone());
            std::fs::create_dir_all(&output).unwrap();
            let db = databases::load_database(std::path::Path::new(&a.db)).unwrap();
            classify::classify(a, &db);
            output
        }
        cli::Commands::Sintax(a) => {
            pool(a.threads);
            let output = a.output_dir.clone().unwrap_or(a.input_dir.clone());
            std::fs::create_dir_all(&output).unwrap();
            let db = databases::load_database(std::path::Path::new(&a.db)).unwrap();
            sintax::sintax(a, &db);
            output
        }
        cli::Commands::Export(a) => {
            merge::export(a);
            a.output_dir.clone()
        }
        _ => panic!("taxonomy: expected classify, sintax or export"),
    };
    let mut names: Vec<String> = std::fs::read_dir(&output_dir)
        .unwrap()
        .map(|e| e.unwrap().file_name().into_string().unwrap())
        .collect();
    names.sort();
    for name in names {
        let text = std::fs::read_to_string(std::path::Path::new(&output_dir).join(&name)).unwrap();
        for line in text.lines() {
            println!("{}	{}", name, line);
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("primitives") => primitives(&args[2..]),
        Some("stats") => stats(),
        Some("markers") => markers(&args[2..]),
        Some("cluster") => cluster(&args[2..]),
        Some("poa") => poa(&args[2..]),
        Some("consensus") => consensus(&args[2..]),
        Some("seedmap") => seedmap(&args[2..]),
        Some("polish") => polish(&args[2..]),
        Some("merge") => merge(&args[2..]),
        Some("chimera") => chimera_stage(&args[2..]),
        Some("em") => em(&args[2..]),
        Some("asv") => asv(&args[2..]),
        Some("taxonomy") => taxonomy_command(&args[2..]),
        _ => {
            eprintln!("usage: savont-upstream-oracle primitives|stats|markers|cluster|poa ...");
            std::process::exit(2);
        }
    }
}
