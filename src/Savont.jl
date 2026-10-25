module Savont

using CodecZlib
using Downloads
using FASTX
using Minimap2
using Printf
using Tar

include("configuration.jl")
include("sequences.jl")
include("io.jl")
include("kmers.jl")
include("statistics.jl")
include("markers.jl")
include("twin_reads.jl")
include("clustering.jl")
include("marker_clustering.jl")
include("poa.jl")
include("consensus.jl")
include("polishing.jl")
include("merging.jl")
include("chimera.jl")
include("quantification.jl")
include("pipeline.jl")
include("taxonomy.jl")
include("classify.jl")
include("sintax.jl")
include("export.jl")
include("download.jl")
include("cli.jl")

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
export POAScoring, POABand, POAGraph, poa_add!, poa_align, poa_consensus
export ConsensusSequence, initial_consensus
export PileupEntry, Pileup, PILEUP_BASE, PILEUP_DELETION, PILEUP_INSERTION, consensus_pileup, top_cluster_indices
export quality_error_rates, analyze_pileup!, is_low_quality, decompress!, polish_consensuses!
export remove_similar_consensuses, adjusted_errors, merge_similar_consensuses
export ChimeraInfo, pairwise_similarities, detect_chimeras, filter_chimeras
export ReadMapping, asv_twin_reads, refine_depths!, per_sample_depths
export write_asv_fasta, write_feature_table, write_final_assignments, savont_asv
export TaxonomyEntry, TaxonomyDatabase, AsvClassification, load_database, database_key, read_fasta_records
export assign_taxonomy, write_species_abundance, write_genus_abundance, write_asv_mappings
export database_mapper, taxon_em, savont_classify, SintaxHit, sintax_kmers, sintax_classify, savont_sintax
export sequence_hash, fuzzy_merge!, savont_export, preprocess_unite, download_database

end
