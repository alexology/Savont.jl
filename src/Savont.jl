module Savont

using CodecZlib
using FASTX

include("configuration.jl")
include("sequences.jl")
include("io.jl")
include("kmers.jl")

export SavontConfig, apply_preset, validate
export ReadRecord, read_records
export reverse_complement, homopolymer_compress, homopolymer_decompress
export homopolymer_compress_with_quality, quality_pool, estimate_sequence_identity
export Blockmer, encode_base, decode_kmer, mm_hash64, rev_hash64
export fmh_seeds, fmh_seeds_positions, minimizer_seeds_positions
export split_kmer_mid, blockmer_kmers, retrieve_masked_kmer, split_kmer

end
