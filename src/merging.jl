# Stage 5: merging of similar consensuses (upstream `alignment.rs`:
# remove_similar_seqs_kmers, calculate_adjusted_errors, merge_similar_consensuses).

const ADAPTER_BUFFER = 25
const ERROR_BUFFER = 35

"""
    remove_similar_consensuses(consensuses; upstream_quirks=false) -> Vector{ConsensusSequence}

Upstream k-mer deduplication: drop a consensus when every minimizer (`w = 10`,
`k = 21`) of its sequence without 25 bases at each end also occurs in one
consensus with more than twice its depth. Consensuses shorter than 100 bases are
dropped as well (upstream skips them when indexing and then discards them).
Order is kept. `upstream_quirks` selects upstream's minimizers (deviation D1) and
compares against zero depth for consensuses already kept, so that, as upstream,
no consensus sorted by decreasing depth is ever removed (deviation D5).
"""
function remove_similar_consensuses(consensuses::AbstractVector{ConsensusSequence}; upstream_quirks::Bool=false)
    index = Dict{UInt64,Vector{Int}}()
    minimizers = Dict{Int,Vector{UInt64}}()
    for (i, c) in enumerate(consensuses)
        length(c.sequence) < 100 && continue
        values, _ = minimizer_seeds_positions(view(c.sequence, (ADAPTER_BUFFER + 1):(length(c.sequence) - ADAPTER_BUFFER)),
                                              10, 21; upstream_quirks)
        for value in values
            push!(get!(Vector{Int}, index, value), i)
        end
        minimizers[i] = values
    end
    # Upstream moves each kept consensus out of the vector (`std::mem::take`),
    # leaving depth 0 behind for the later comparisons (deviation D5).
    depths = [c.depth for c in consensuses]
    kept = ConsensusSequence[]
    for i in eachindex(consensuses)
        values = get(minimizers, i, nothing)
        values === nothing && continue
        deeper = Set{Int}()
        for (k, value) in enumerate(values)
            ids = index[value]
            if k == 1
                for id in ids
                    depths[id] ÷ 2 > depths[i] && push!(deeper, id)
                end
            else
                intersect!(deeper, ids)
            end
        end
        if isempty(deeper)
            push!(kept, consensuses[i])
            upstream_quirks && (depths[i] = 0)
        end
    end
    return kept
end

# Upstream `has_homopolymer_context` with 0-based `position`: a run of three equal
# bases within `window` bases of the position.
function _homopolymer_context(sequence, position::Int, window::Int)
    isempty(sequence) && return false
    first = max(position - window, 0)
    last = min(position + window + 1, length(sequence))
    last <= first + 2 && return false
    for i in first:(last - 3)
        if i + 2 < length(sequence) && sequence[i + 1] == sequence[i + 2] == sequence[i + 3]
            return true
        end
    end
    return false
end

"""
    adjusted_errors(cigar, query, target, query_start, target_start) -> Int

Upstream `calculate_adjusted_errors` (0-based starts, minimap2 CIGAR words):
mismatches other than with `N`, and indels not within two bases of a homopolymer
run of three (one error below 10 bases, the length otherwise), both counted only
more than 35 bases from the query (mismatches, insertions) or target (deletions)
ends.
"""
function adjusted_errors(cigar::AbstractVector{UInt32}, query, target, query_start::Int, target_start::Int)
    errors = 0
    q = query_start; t = target_start
    for word in cigar
        n = Int(word >> 4); op = word & 0xf
        if op == 0
            for _ in 1:n
                if q < length(query) && t < length(target)
                    a = query[q + 1]; b = target[t + 1]
                    if a != b && a != UInt8('N') && b != UInt8('N') && q > ERROR_BUFFER && q + ERROR_BUFFER < length(query)
                        errors += 1
                    end
                end
                q += 1; t += 1
            end
        elseif op == 1
            if !(_homopolymer_context(query, q, 2) || _homopolymer_context(target, t, 2)) &&
               q > ERROR_BUFFER && q + n + ERROR_BUFFER < length(query)
                errors += n < 10 ? 1 : n
            end
            q += n
        elseif op == 2
            if !(_homopolymer_context(query, q, 2) || _homopolymer_context(target, t, 2)) &&
               t > ERROR_BUFFER && t + n + ERROR_BUFFER < length(target)
                errors += n < 10 ? 1 : n
            end
            t += n
        else
            @warn "Unexpected CIGAR operation: $op"
        end
    end
    return errors
end

# minimap2 as configured by upstream for stage 5 through the minimap2 crate:
# `lrhq().with_cigar().with_index(fasta)` (mapping options updated from the index),
# then `set_no_diag()` and `best_n = 75`.
function _consensus_mapper(sequences::AbstractVector{<:AbstractVector{UInt8}})
    index_options, map_options = Minimap2.preset("lr:hq")
    index = Minimap2.MinimapIndex(sequences, index_options; names=["c$i" for i in eachindex(sequences)])
    Minimap2.update_options!(map_options, index)
    map_options.flag |= Minimap2.F_CIGAR | Minimap2.F_NO_DIAG
    map_options.best_n = 75
    return index, map_options
end

"""
    merge_similar_consensuses(consensuses, config) -> Vector{ConsensusSequence}

Stage 5. After `remove_similar_consensuses`, every consensus (decompressed
sequence) is mapped to all the others with minimap2 (`lr:hq`). A mapping covering
at least 3/4 of the query with NM ≤ 30 gives `e` = min(NM, `adjusted_errors`).
The query merges into a target when their depth ratio (either way) is below
0.5^(0.75e + 1.25), or 0.999999 when e = 0; identical consensuses of equal depth
merge into the earlier one. A query merges into its deepest valid target (then
fewest errors, then index); shallower valid targets merge into the query unless
already assigned; chains are followed to their end. Merged consensuses get the
union of the clusters (sorted) as their cluster and its size as depth, lose
their low-quality positions, and are sorted by decreasing depth, then sequence,
then index. Low-quality consensuses are not needed: upstream maps them only to
increase a depth field that this function resets. If no consensus survives
`remove_similar_consensuses` the result is empty (upstream panics; deviation D7).
"""
function merge_similar_consensuses(polished::AbstractVector{ConsensusSequence}, config::SavontConfig)
    isempty(polished) && return ConsensusSequence[]
    consensuses = remove_similar_consensuses(polished; upstream_quirks=config.upstream_quirks)
    # Deviation D7: every consensus can be shorter than 100 bases
    if isempty(consensuses)
        config.upstream_quirks &&
            error("Failed to create aligner: Index File is empty (upstream panics here, deviation D7)")
        return ConsensusSequence[]
    end
    for c in consensuses
        c.decompressed === nothing && decompress!(c)
    end
    sequences = [c.decompressed for c in consensuses]
    index, options = _consensus_mapper(sequences)
    # (target, adjusted errors) per query, in region order.
    mappings = [Tuple{Int,Int}[] for _ in consensuses]
    Threads.@threads :dynamic for query in eachindex(consensuses)
        sequence = sequences[query]
        regions, _ = Minimap2.map_query(index, sequence, options)
        for r in regions
            r.extra === nothing && continue
            target = Int(r.rid) + 1
            nm = Int(r.blen) - Int(r.mlen) + Int(r.extra.n_ambi)
            (r.qe - r.qs < length(sequence) * 3 ÷ 4 || nm > 30) && continue
            errors = if isempty(r.extra.cigar)
                nm
            elseif r.rev
                adjusted_errors(r.extra.cigar, reverse_complement(sequence), sequences[target], length(sequence) - Int(r.qe), Int(r.rs))
            else
                adjusted_errors(r.extra.cigar, sequence, sequences[target], Int(r.qs), Int(r.rs))
            end
            push!(mappings[query], (target, min(errors, nm)))
        end
    end
    merge_map = Dict{Int,Int}()
    depth(i) = consensuses[i].depth
    for query in eachindex(consensuses)
        to_deeper = Tuple{Int,Int,Int}[]       # (target depth, errors, target)
        for (target, errors) in mappings[query]
            query == target && continue
            relative = depth(query) / depth(target)
            threshold = 0.5^(errors * 0.75 + 1.25)
            valid = if errors == 0 && depth(query) == depth(target)
                query > target
            else
                errors == 0 && (threshold = 0.999999)
                relative < threshold || 1.0 / relative < threshold
            end
            valid || continue
            if depth(target) == depth(query)
                errors == 0 && query > target && (merge_map[query] = target)
            elseif depth(target) > depth(query)
                push!(to_deeper, (depth(target), errors, target))
            else
                haskey(merge_map, target) || (merge_map[target] = query)
            end
        end
        if !isempty(to_deeper)
            sort!(to_deeper; by=x -> (-x[1], x[2], x[3]))
            merge_map[query] = to_deeper[1][3]
        end
    end
    clusters = [copy(c.cluster) for c in consensuses]
    for query in eachindex(consensuses)
        haskey(merge_map, query) || continue
        target = merge_map[query]
        while haskey(merge_map, target)
            target = merge_map[target]
        end
        append!(clusters[target], clusters[query]); empty!(clusters[query])
    end
    merged = ConsensusSequence[]
    for (c, cluster) in zip(consensuses, clusters)
        isempty(cluster) && continue
        sort!(cluster)
        push!(merged, decompress!(ConsensusSequence(c.id, c.sequence, c.hp_lengths, length(cluster), cluster)))
    end
    return sort!(merged; by=c -> (-c.depth, c.sequence, c.id))
end
