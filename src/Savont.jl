module Savont

using CodecZlib
using FASTX

include("configuration.jl")
include("sequences.jl")
include("io.jl")
include("kmers.jl")
include("statistics.jl")
include("markers.jl")
include("twin_reads.jl")
include("clustering.jl")
include("marker_clustering.jl")

export SavontConfig, apply_preset, validate
export ReadRecord, read_records, foreach_record
export reverse_complement, homopolymer_compress, homopolymer_decompress
export homopolymer_compress_with_quality, quality_pool, estimate_sequence_identity
export Blockmer, encode_base, decode_kmer, mm_hash64, rev_hash64
export fmh_seeds, fmh_seeds_positions, minimizer_seeds_positions
export split_kmer_mid, blockmer_kmers, retrieve_masked_kmer, split_kmer
export binomial_test, fishers_exact
export KmerCount, SnpmerInfo, BlockmerInfo, KmerGlobalInfo, snpmer_kmers
export count_split_kmers, find_snpmers, find_blockmers, discover_markers
export TwinRead, twin_read, lsh_signatures, build_twin_reads, cluster_reads_by_kmers
export snpmer_markers, blockmer_markers, ConsensusMarker, greedy_snpmer_clusters, cluster_reads_by_snpmers

end
