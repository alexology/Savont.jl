//! Differential oracle built from the pinned upstream savont 0.8.0 sources.
//! Source commit: 5a5a1ce9053a867abc09e991f61b8ec7fb4fc6ef.
//!
//! The upstream modules are compiled unchanged; this file only drives them and
//! prints tab-separated records for test/rust_oracle.jl. Driver code copied from
//! main.rs (which cannot be included) is marked as such.
#![allow(dead_code, unused_imports, unused_variables, unused_mut)]

#[path = "../../../reference/savont/src/asv_cluster.rs"]
mod asv_cluster;
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

/// cli.rs only needs the download keywords; the full module pulls in taxonomy.
mod databases {
    pub const KEYWORDS: &[&str] = &["emu-1", "silva-138.2", "greengenes2-2024.09", "unite-10.0"];
}

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
    println!("low_polymorphism\t{}", u8::from(fraction_without > 0.75));
    // As main.rs: auto-enable low-polymorphism mode.
    if fraction_without > 0.75 && !cluster_args.low_polymorphism {
        cluster_args.low_polymorphism = true;
    }

    for (index, read) in twin_reads.iter().enumerate() {
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
    for members in clusters.iter() {
        println!("cluster\t{}", csv(members));
    }
    let final_clusters =
        asv_cluster::cluster_reads_by_snpmers(&twin_reads, &clusters, &cluster_args, &directory);
    // Greedy SNPmer clusters before reclustering, as written by upstream.
    let greedy = directory.join("snpmer_clusters_before_reclust2.5.tsv");
    if let Ok(text) = std::fs::read_to_string(&greedy) {
        for line in text.lines().skip(1) {
            let fields = line.split('\t').collect::<Vec<_>>();
            println!("greedy\t{}\t{}", fields[0], fields[4]);
        }
    }
    std::fs::remove_dir_all(&directory).unwrap();
    for members in final_clusters.iter() {
        println!("final\t{}", csv(members));
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("primitives") => primitives(&args[2..]),
        Some("stats") => stats(),
        Some("markers") => markers(&args[2..]),
        Some("cluster") => cluster(&args[2..]),
        _ => {
            eprintln!("usage: savont-upstream-oracle primitives|stats|markers|cluster ...");
            std::process::exit(2);
        }
    }
}
