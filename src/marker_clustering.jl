# Phase 3C, second part: polymorphic-marker clustering (upstream
# `asv_cluster::cluster_reads_by_snpmers` and `recluster_using_consensus_reps`).

# Upstream switches from the SNPmer index to a scan over representatives above
# this many representatives in one k-mer cluster.
const ITERATIVE_REPRESENTATIVES = 1000

function _encode_kmer(sequence::AbstractVector{UInt8}, position::Integer, length::Integer)
    kmer = UInt64(0)
    for i in (position + 1):(position + length)
        kmer = (kmer << 2) | encode_base(sequence[i])
    end
    return kmer
end
# The same from a read's 2-bit sequence.
function _encode_kmer(read::TwinRead, position::Integer, length::Integer)
    packed = getfield(read, :packed_sequence)
    kmer = UInt64(0)
    for i in (position + 1):(position + length)
        kmer = (kmer << 2) | _packed_code(packed, i)
    end
    return kmer
end

function _reverse_complement_kmer(kmer::UInt64, length::Integer)
    result = UInt64(0)
    for _ in 1:length
        result = (result << 2) | (UInt64(3) - (kmer & UInt64(3)))
        kmer >>= 2
    end
    return result
end

"""
    snpmer_markers(read::TwinRead) -> Vector{Tuple{UInt32,UInt64}}

`(position, k-mer)` for the retained SNPmer positions, with the k-mer re-read from
`read.sequence` (upstream `TwinRead::snpmers_vec`): the reverse complement is used
when its middle-masked value is smaller, otherwise the forward k-mer.
"""
function snpmer_markers(read::TwinRead)
    k = Int(read.k); mask = ~(UInt64(3) << (k - 1))
    return map(read.snpmer_positions) do position
        forward = _encode_kmer(read, position, k)
        reverse = _reverse_complement_kmer(forward, k)
        (position, (reverse & mask) < (forward & mask) ? reverse : forward)
    end
end

"""
    blockmer_markers(read::TwinRead) -> Vector{Tuple{UInt32,UInt64}}

`(position, blockmer)` re-read from `read.sequence` in the orientation that
matched (upstream `TwinRead::blockmers_vec`).
"""
function blockmer_markers(read::TwinRead)
    length_kl = Int(read.k) + Int(read.l)
    return map(zip(read.blockmer_positions, read.blockmer_canonical)) do (position, forward)
        kmer = _encode_kmer(read, position, length_kl)
        (position, forward ? kmer : _reverse_complement_kmer(kmer, length_kl))
    end
end

# Blockmer matches and mismatches between two reads; for repeated anchors in
# `other`, the last occurrence counts (upstream `compare_blockmers`).
function _compare_blockmers(read::TwinRead, other::TwinRead, l::Integer)
    anchors = Dict{UInt64,UInt64}()
    for (_, kmer) in blockmer_markers(other)
        anchors[kmer >> (2l)] = kmer
    end
    matches = 0; mismatches = 0
    for (_, kmer) in blockmer_markers(read)
        other_kmer = get(anchors, kmer >> (2l), nothing)
        other_kmer === nothing && continue
        other_kmer == kmer ? (matches += 1) : (mismatches += 1)
    end
    return matches, mismatches
end

function _best_representative_indexed(read_index, reads, index, representative_size, k, l, use_blockmers)
    mask = ~(UInt64(3) << (k - 1))
    statistics = Dict{Int,Tuple{Int,Int}}()
    for kmer in reads[read_index].snpmer_kmers
        for (candidate, candidate_kmer) in get(index, kmer & mask, ())
            matches, mismatches = get(statistics, candidate, (0, 0))
            statistics[candidate] = candidate_kmer == kmer ? (matches + 1, mismatches) : (matches, mismatches + 1)
        end
    end
    candidates = [(-matches, representative_size[candidate], candidate)
                  for (candidate, (matches, mismatches)) in statistics if mismatches == 0 && matches > 0]
    isempty(candidates) && return 0
    sort!(candidates)
    use_blockmers || return candidates[1][3]
    checked = [(candidate, _compare_blockmers(reads[read_index], reads[candidate], l)...) for (_, _, candidate) in candidates]
    sort!(checked; by=x -> (x[3], -x[2], x[1]))
    return checked[1][3] > 1 ? 0 : checked[1][1]
end

function _best_representative_iterative(read_index, reads, representatives, k, l, use_blockmers)
    mask = ~(UInt64(3) << (k - 1))
    query = Dict(kmer & mask => kmer for kmer in reads[read_index].snpmer_kmers)
    for candidate in representatives
        matches = 0; mismatches = 0
        for kmer in reads[candidate].snpmer_kmers
            query_kmer = get(query, kmer & mask, nothing)
            query_kmer === nothing && continue
            query_kmer == kmer ? (matches += 1) : (mismatches += 1)
        end
        mismatches == 0 && matches > 0 || continue
        use_blockmers && _compare_blockmers(reads[read_index], reads[candidate], l)[2] > l && continue
        return candidate
    end
    return 0
end

_size_order(cluster) = (-length(cluster), first(cluster))

# Greedy SNPmer clustering of one k-mer cluster (members in increasing order).
function _greedy_snpmer_clusters(kmer_cluster, reads, config)
    k = config.kmer_size; l = config.blockmer_length; mask = ~(UInt64(3) << (k - 1))
    index = Dict{UInt64,Vector{Tuple{Int,UInt64}}}()
    representatives = Int[]
    representative_size = Dict{Int,Int}()
    assignment = Dict{Int,Int}()
    for read_index in kmer_cluster
        best = if length(representatives) > ITERATIVE_REPRESENTATIVES
            _best_representative_iterative(read_index, reads, representatives, k, l, config.use_blockmers)
        else
            _best_representative_indexed(read_index, reads, index, representative_size, k, l, config.use_blockmers)
        end
        if best == 0
            push!(representatives, read_index)
            for kmer in reads[read_index].snpmer_kmers
                push!(get!(Vector{Tuple{Int,UInt64}}, index, kmer & mask), (read_index, kmer))
            end
            representative_size[read_index] = 1
            best = read_index
        else
            representative_size[best] += 1
        end
        assignment[read_index] = best
    end
    members = Dict{Int,Vector{Int}}()
    for (read_index, representative) in assignment
        push!(get!(Vector{Int}, members, representative), read_index)
    end
    clusters = sort!([sort!(cluster) for cluster in values(members)]; by=_size_order)
    return filter!(cluster -> length(cluster) >= config.min_cluster_size, clusters)
end

"""
    ConsensusMarker

Consensus of one polymorphic marker in a cluster: median position, the context
(`splitmer`: middle-masked k-mer for SNPmers, anchor k-mer for blockmers), the most
frequent full k-mer and its count.
"""
struct ConsensusMarker
    position::UInt32
    splitmer::UInt64
    kmer::UInt64
    count::Int
end

function _consensus_markers(cluster, reads, use_blockmers::Bool, k, l)
    mask = ~(UInt64(3) << (k - 1))
    data = Dict{UInt64,Dict{UInt64,Tuple{Int,Vector{UInt32}}}}()
    for read_index in cluster
        markers = use_blockmers ? blockmer_markers(reads[read_index]) : snpmer_markers(reads[read_index])
        for (position, kmer) in markers
            splitmer = use_blockmers ? kmer >> (2l) : kmer & mask
            alleles = get!(Dict{UInt64,Tuple{Int,Vector{UInt32}}}, data, splitmer)
            count, positions = get!(() -> (0, UInt32[]), alleles, kmer)
            push!(positions, position)
            alleles[kmer] = (count + 1, positions)
        end
    end
    minimum_count = max(length(cluster) ÷ 6, 1)
    consensus = ConsensusMarker[]
    for (splitmer, alleles) in data
        # Most frequent allele; ties go to the smallest k-mer.
        kmer, (count, positions) = first(sort!(collect(alleles); by=x -> (-x[2][1], x[1])))
        count >= minimum_count || continue
        push!(consensus, ConsensusMarker(sort(positions)[length(positions) ÷ 2 + 1], splitmer, kmer, count))
    end
    return sort!(consensus; by=x -> (x.position, x.splitmer))
end

function _compare_consensus(first_consensus, second_consensus)
    second = Dict(x.splitmer => x.kmer for x in second_consensus)
    matches = 0; mismatches = 0
    for marker in first_consensus
        kmer = get(second, marker.splitmer, nothing)
        kmer === nothing && continue
        kmer == marker.kmer ? (matches += 1) : (mismatches += 1)
    end
    return matches, mismatches
end

# Upstream `are_consensus_concordant`, including its threshold
# min(length(first), max(length(second), 2)).
function _concordant(first_consensus, second_consensus)
    matches, mismatches = _compare_consensus(first_consensus, second_consensus)
    return mismatches == 0 && matches >= min(length(first_consensus), max(length(second_consensus), 2))
end

# One merge round within a k-mer group (upstream `recluster_one_round`). Clusters
# are visited by decreasing size; each absorbs every later cluster whose consensus
# is concordant with its own consensus as computed before the merges.
function _merge_round(clusters, reads, use_blockmers, k, l)
    entries = [(copy(cluster), _consensus_markers(cluster, reads, use_blockmers, k, l))
               for cluster in clusters if !isempty(cluster)]
    sort!(entries; by=x -> _size_order(x[1]))
    merged = falses(length(entries))
    result = Vector{Vector{Int}}()
    merges = 0
    for i in eachindex(entries)
        merged[i] && continue
        consensus_i = entries[i][2]
        for j in (i + 1):length(entries)
            merged[j] && continue
            consensus_j = entries[j][2]
            concordant = _concordant(consensus_i, consensus_j) && _concordant(consensus_j, consensus_i)
            matches, mismatches = _compare_consensus(consensus_i, consensus_j)
            larger = max(length(entries[i][1]), length(entries[j][1]))
            smaller = min(length(entries[i][1]), length(entries[j][1]))
            shortest = min(length(consensus_i), length(consensus_j))
            if mismatches == 0 && matches > shortest * 0.975 && larger ÷ smaller > 50
                concordant = true
            end
            if mismatches == 0 && larger ÷ smaller > 500 && smaller <= 2
                concordant = true
            end
            if concordant
                append!(entries[i][1], entries[j][1])
                merged[j] = true
                merges += 1
            end
        end
        push!(result, sort!(entries[i][1]))
    end
    return sort!(result; by=_size_order), merges
end

# Reassign every read to the cluster whose consensus it matches best (fewest
# mismatches, then most matches, then earliest cluster); upstream
# `reassign_reads_to_best_cluster`.
function _reassign_round(clusters, reads, use_blockmers, k, l, min_cluster_size)
    mask = ~(UInt64(3) << (k - 1))
    lookups = [Dict(x.splitmer => x.kmer for x in _consensus_markers(cluster, reads, use_blockmers, k, l))
               for cluster in clusters]
    result = [Int[] for _ in clusters]
    reassignments = 0
    for (cluster_index, cluster) in enumerate(clusters), read_index in cluster
        markers = use_blockmers ? blockmer_markers(reads[read_index]) : snpmer_markers(reads[read_index])
        best = cluster_index; best_mismatches = typemax(Int); best_matches = 0
        for (candidate, lookup) in enumerate(lookups)
            matches = 0; mismatches = 0
            for (_, kmer) in markers
                consensus_kmer = get(lookup, use_blockmers ? kmer >> (2l) : kmer & mask, nothing)
                consensus_kmer === nothing && continue
                consensus_kmer == kmer ? (matches += 1) : (mismatches += 1)
            end
            if mismatches < best_mismatches || (mismatches == best_mismatches && matches > best_matches)
                best = candidate; best_mismatches = mismatches; best_matches = matches
            end
        end
        push!(result[best], read_index)
        best == cluster_index || (reassignments += 1)
    end
    filter!(cluster -> !isempty(cluster) && length(cluster) >= min_cluster_size, result)
    foreach(sort!, result)
    return result, reassignments
end

"""
    greedy_snpmer_clusters(reads, kmer_clusters, config) -> Vector{Vector{Vector{Int}}}

The greedy pass of [`cluster_reads_by_snpmers`](@ref) alone, one vector of clusters
per primary cluster (upstream `snpmer_clusters_before_reclust2.5.tsv`).
"""
greedy_snpmer_clusters(reads::AbstractVector{TwinRead}, kmer_clusters, config::SavontConfig) =
    _tmap(cluster -> _greedy_snpmer_clusters(cluster, reads, config), kmer_clusters)

"""
    cluster_reads_by_snpmers(reads, kmer_clusters, config; low_polymorphism=config.low_polymorphism)
        -> Vector{Vector{Int}}

Polymorphic-marker clustering within each primary cluster (upstream stage 3).

1. Greedy pass: reads are visited in order; a read joins the representative whose
   SNPmers (`snpmer_kmers`) share at least one allele and contradict none, preferring
   more shared alleles and then the smaller cluster. With `use_blockmers`, the chosen
   representative must also have at most one blockmer mismatch. Otherwise the read
   becomes a representative. Reads without SNPmers therefore always start their own
   cluster, and clusters smaller than `min_cluster_size` are dropped.
2. Reclustering, at most `max_iterations_recluster` rounds within each primary
   cluster, until a round makes no merge: clusters whose consensus markers (alleles
   seen in at least 1/6 of the reads) agree are merged, with additional merges for
   very unequal sizes; then every read moves to the cluster whose consensus it
   matches best, and clusters below `min_cluster_size` are dropped. Consensus
   markers are SNPmers, or blockmers when `use_blockmers`.

With `low_polymorphism` the primary clusters are returned unchanged (sizes at least
`min_cluster_size`). Otherwise, unless `config.upstream_quirks`, the same applies to
each primary cluster in which more than 75% of the reads have no SNPmer: upstream
applies this rule only to the whole data set, so a low-polymorphism amplicon among
polymorphic ones would be discarded (deviation D3 in `docs/deviations.md`).
Clusters are sorted vectors of 1-based read indices, ordered by decreasing size and
then by first member.
"""
function cluster_reads_by_snpmers(reads::AbstractVector{TwinRead}, kmer_clusters::AbstractVector{<:AbstractVector{<:Integer}},
                                  config::SavontConfig; low_polymorphism::Bool=config.low_polymorphism)
    if low_polymorphism
        kept = [collect(Int, cluster) for cluster in kmer_clusters if length(cluster) >= config.min_cluster_size]
        return sort!(kept; by=_size_order)
    end
    k = config.kmer_size; l = config.blockmer_length; use_blockmers = config.use_blockmers
    passed = Vector{Vector{Int}}()
    if !config.upstream_quirks
        polymorphic = Vector{Vector{Int}}()
        for cluster in kmer_clusters
            without = count(i -> isempty(reads[i].snpmer_positions), cluster)
            push!(without / length(cluster) > 0.75 ? passed : polymorphic, collect(Int, cluster))
        end
        kmer_clusters = polymorphic
    end
    groups = greedy_snpmer_clusters(reads, kmer_clusters, config)
    for _ in 1:config.max_iterations_recluster
        # Groups (primary clusters) are independent.
        rounds = _tmap(groups) do group
            merged, merges = _merge_round(group, reads, use_blockmers, k, l)
            first(_reassign_round(merged, reads, use_blockmers, k, l, config.min_cluster_size)), merges
        end
        groups = first.(rounds)
        sum(last, rounds; init=0) == 0 && break
    end
    final = sort!(vcat([cluster for group in groups for cluster in group], passed); by=_size_order)
    return filter!(cluster -> length(cluster) >= config.min_cluster_size, final)
end
