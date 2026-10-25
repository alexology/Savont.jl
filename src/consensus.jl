# Phase 3D: initial consensus of each cluster (upstream
# `alignment::align_and_consensus`, stage 4a).

const MAX_SEQS_CONSENSUS_POA = 75

"""
    ConsensusSequence

Consensus of one cluster: the (optionally homopolymer-compressed) `sequence`, the
homopolymer lengths of its bases (1 until pileups estimate them), the number of
reads in the cluster (`depth`), the 1-based index of the cluster (`id`) and its
read indices. Polishing (stage 4b) adds the 1-based `low_quality_positions` and
the `decompressed` sequence without masked ends; chimera detection (stage 6) the
`chimera_score`; quantification (stage 7) the counts of reads assigned to it
alone, together with other ASVs, and with NM ≤ 10, and with pooled samples the
`per_sample_depths`.
"""
mutable struct ConsensusSequence
    id::Int
    sequence::Vector{UInt8}
    hp_lengths::Vector{UInt8}
    depth::Int
    cluster::Vector{Int}
    low_quality_positions::Vector{Int}
    decompressed::Union{Nothing,Vector{UInt8}}
    chimera_score::Union{Nothing,Int}
    unambiguous_count::Union{Nothing,Int}
    ambiguous_count::Union{Nothing,Int}
    close_count::Union{Nothing,Int}
    per_sample_depths::Vector{Int}
end
ConsensusSequence(id, sequence, hp_lengths, depth, cluster) =
    ConsensusSequence(id, sequence, hp_lengths, depth, cluster, Int[], nothing, nothing, nothing, nothing, nothing, Int[])

# minimap2 as configured by upstream through the minimap2 Rust crate:
# `Aligner::builder().map_ont().with_cigar().with_seq(seed)`, i.e. the map-ont
# preset (or `preset`), an index of the single sequence named "N/A", mid_occ fixed
# at 1000 (no index-based update) and queries mapped without a name.
# `bucket_bits` is the index option upstream lowers for one-off references.
function _seed_mapper(seed::AbstractVector{UInt8}; preset::String="map-ont", bucket_bits::Int=14)
    index_options, map_options = Minimap2.preset(preset)
    index_options.bucket_bits = bucket_bits
    index = Minimap2.MinimapIndex([seed], index_options; names=["N/A"])
    map_options.mid_occ = 1000
    map_options.flag |= Minimap2.F_CIGAR
    return index, map_options
end

# Qualities as upstream rebuilds them from the 4-base bins: Phred+33 bytes 3 * code
# + 33 repeated over each bin and cut or padded (with the last value) to the read
# length; 33 for reads without qualities.
function _expanded_qualities(read::TwinRead)
    n = read.base_length
    binned = read.binned_qualities === nothing ? fill(UInt8(33), n) : UInt8[3 * c + 33 for c in read.binned_qualities]
    expanded = repeat(binned; inner=QUALITY_SEQ_BIN)
    if length(expanded) > n
        resize!(expanded, n)
    elseif length(expanded) < n
        append!(expanded, fill(expanded[end], n - length(expanded)))
    end
    return binned, expanded
end

"""
    initial_consensus(reads, clusters, config) -> Vector{ConsensusSequence}

Partial-order-alignment consensus of each cluster (upstream stage 4a). The seed is
the read at the 90th percentile of length; up to 75 reads with the highest mean
binned base accuracy are mapped to it with minimap2 (`map-ont`) and, if they map,
oriented on the strand of their best mapping. The seed and those reads, in cluster
order, are aligned with spoa's overlap alignment (match 3, mismatch −8, gap −6;
abPOA band of 10% of the length plus the largest length deviation, unless
`no_band`), weighted by their qualities. With `use_hpc`, reads and consensus are
homopolymer-compressed. Consensuses shorter than 40 bases are dropped. The result
is sorted by decreasing depth, then cluster index. Clusters are processed in
parallel on Julia's threads; the result does not depend on their number.
"""
function initial_consensus(reads::AbstractVector{TwinRead}, clusters::AbstractVector{<:AbstractVector{<:Integer}},
                           config::SavontConfig)
    results = Vector{Union{Nothing,ConsensusSequence}}(nothing, length(clusters))
    # Clusters are independent; each task keeps its own minimap2 and POA buffers.
    Threads.@threads :dynamic for cluster_index in eachindex(clusters)
        results[cluster_index] = _cluster_consensus(reads, clusters[cluster_index], cluster_index, config)
    end
    consensuses = ConsensusSequence[c for c in results if c !== nothing]
    return sort!(consensuses; by=c -> (-c.depth, c.id))
end

function _cluster_consensus(reads::AbstractVector{TwinRead}, cluster::AbstractVector{<:Integer}, cluster_index::Int,
                            config::SavontConfig)
    isempty(cluster) && return nothing
    scoring = POAScoring(3, -8, -6)
    sequences = Vector{UInt8}[]; qualities = Vector{UInt8}[]; mean_accuracy = Float64[]
    for read_index in cluster
        read = reads[read_index]
        binned, expanded = _expanded_qualities(read)
        total = 0.0
        for q in binned
            total += 1.0 - 10.0^(-(Float64(q) - 33) / 10.0)
        end
        push!(mean_accuracy, isempty(binned) ? 1.0 : total / length(binned))
        push!(sequences, read.sequence); push!(qualities, expanded)
    end
    n = length(sequences)
    by_length = sortperm(length.(sequences))                  # stable, as sort_by_key
    seed = by_length[unsafe_trunc(Int, n * 0.9) + 1]
    by_accuracy = sortperm(mean_accuracy; rev=true)           # stable, as sort_by
    index, map_options = _seed_mapper(sequences[seed])
    mapped = Tuple{Int,Vector{Minimap2.Region}}[]
    for i in by_accuracy[1:min(MAX_SEQS_CONSENSUS_POA, n)]
        i == seed && continue
        regions, _ = Minimap2.map_query(index, sequences[i], map_options)
        push!(mapped, (i, regions))
    end
    sort!(mapped; by=first)
    aligned = [sequences[seed]]; aligned_qualities = [qualities[seed]]
    for (i, regions) in mapped
        (i == seed || isempty(regions)) && continue
        if regions[1].rev
            push!(aligned, reverse_complement(sequences[i])); push!(aligned_qualities, reverse(qualities[i]))
        else
            push!(aligned, sequences[i]); push!(aligned_qualities, qualities[i])
        end
        length(aligned) > MAX_SEQS_CONSENSUS_POA && break
    end
    poa_sequences = Vector{UInt8}[]; poa_weights = Vector{UInt8}[]
    for (sequence, quality) in zip(aligned, aligned_qualities)
        compressed, compressed_quality, _ = homopolymer_compress_with_quality(sequence, quality, config.use_hpc)
        push!(poa_sequences, compressed); push!(poa_weights, compressed_quality)
    end
    mean_length = sum(length, poa_sequences) ÷ length(poa_sequences)
    deviation = maximum(abs(mean_length - length(s)) for s in poa_sequences)
    band = config.no_band ? nothing : POABand(deviation, 0.1f0)
    graph = POAGraph()
    for (sequence, weights) in zip(poa_sequences, poa_weights)
        alignment, _ = poa_align(graph, sequence, scoring; band, isa=config.poa_isa)
        poa_add!(graph, alignment, sequence, weights)
    end
    consensus, _ = homopolymer_compress(poa_consensus(graph), config.use_hpc; upstream_quirks=config.upstream_quirks)
    if length(consensus) < 40
        @warn "Consensus of cluster $cluster_index is too short ($(length(consensus)) bases); skipped"
        return nothing
    end
    return ConsensusSequence(cluster_index, consensus, ones(UInt8, length(consensus)), length(cluster), collect(Int, cluster))
end
