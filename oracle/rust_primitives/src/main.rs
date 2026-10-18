mod upstream_primitives;

use upstream_primitives::{
    blockmer_kmers, fmh_seeds, fmh_seeds_positions, minimizer_seeds_positions,
    mm_hash64, rev_hash_64, split_kmer_mid,
    estimate_sequence_identity, homopolymer_compress, homopolymer_compress_with_quality,
    quality_pool, reverse_complement,
};

fn csv(values: &[u64]) -> String {
    values.iter().map(u64::to_string).collect::<Vec<_>>().join(",")
}

fn main() {
    for value in [0_u64, 1, 3, 0x123456789abcdef0, u64::MAX] {
        let hashed = mm_hash64(value);
        println!("hash\t{value}\t{hashed}\t{}", rev_hash_64(hashed));
    }

    let sequence = b"ACGTACGTACGT";
    let mut hashes = Vec::new();
    fmh_seeds(sequence, &mut hashes, 1, 5);
    println!("fmh\t{}", csv(&hashes));

    let mut kmers = Vec::new();
    let mut positions = Vec::new();
    fmh_seeds_positions(sequence, &mut kmers, &mut positions, 1, 5);
    println!("fmhp\t{}\t{}", csv(&kmers), csv(&positions));

    let split = split_kmer_mid(sequence.to_vec(), Some(vec![33; sequence.len()]), 5, 25);
    println!("split\t{}", csv(&split));

    let blocks = blockmer_kmers(b"AACGTACGTT".to_vec(), None, 5, 2, 25);
    println!(
        "blocks\t{}",
        blocks.iter().map(|x| format!("{}:{}", x.kmer, u8::from(x.is_forward))).collect::<Vec<_>>().join(",")
    );

    let mut minimizers = Vec::new();
    let mut minimizer_positions = Vec::new();
    minimizer_seeds_positions(sequence, &mut minimizers, &mut minimizer_positions, 3, 5);
    println!("minimizer\t{}\t{}", csv(&minimizers), csv(&minimizer_positions));
    println!("revcomp\t{}", String::from_utf8(reverse_complement(b"AaCGTn?")).unwrap());

    let (hpc, lengths) = homopolymer_compress(b"AAACCG", true);
    println!("hpc\t{}\t{}", String::from_utf8(hpc).unwrap(), lengths.iter().map(u8::to_string).collect::<Vec<_>>().join(","));
    let (_, qualities, _) = homopolymer_compress_with_quality(b"AAAC", &[70, 60, 65, 50], true);
    println!("hpcq\t{}", qualities.iter().map(u8::to_string).collect::<Vec<_>>().join(","));
    println!("pool\t{}", quality_pool(vec![9, 8, 7, 6, 5, 4, 3]).iter().map(u8::to_string).collect::<Vec<_>>().join(","));
    println!("identity\t{:.17}", estimate_sequence_identity(&[73]).unwrap());
}
