# Stage 6: chimera detection (upstream `chimera.rs`).

const PCR_SLACK = 15

# Rust `f64 as i64`: truncation, saturating at the bounds, NaN to 0.
_as_i64(x::Float64) = isnan(x) ? 0 : x <= -0x1p63 ? typemin(Int) : x >= 0x1p63 ? typemax(Int) : trunc(Int, x)

"""
    ChimeraInfo

A chimera call: indices (into the consensus vector) of the query and of its left
and right parents (the same for a single-parent call), the matched lengths, the
query length and the fraction of the query covered by the two matches.
"""
struct ChimeraInfo
    query::Int
    left_parent::Int
    right_parent::Int
    left_match_length::Int
    right_match_length::Int
    query_length::Int
    coverage_fraction::Float64
end

"""
    pairwise_similarities(consensuses) -> Dict{Tuple{Int,Int},Float64}

Identity 1 − NM / query span of the first `lr:hq` mapping of consensus `j` to
consensus `i`, for i < j unless `i` is more than 25 times deeper than `j`.
Keyed `(j, i)`, as upstream.
"""
function pairwise_similarities(consensuses::AbstractVector{ConsensusSequence})
    local_results = [Tuple{Tuple{Int,Int},Float64}[] for _ in consensuses]
    Threads.@threads :dynamic for i in eachindex(consensuses)
        index, options = _seed_mapper(consensuses[i].decompressed; preset="lr:hq")
        for j in (i + 1):length(consensuses)
            consensuses[i].depth > consensuses[j].depth * 25 && continue
            regions, _ = Minimap2.map_query(index, consensuses[j].decompressed, options)
            isempty(regions) && continue
            r = regions[1]
            r.extra === nothing && continue
            span = Float64(r.qe - r.qs)
            nm = Int(r.blen) - Int(r.mlen) + Int(r.extra.n_ambi)
            push!(local_results[i], ((j, i), span > 0 ? 1.0 - nm / span : 0.0))
        end
    end
    return Dict{Tuple{Int,Int},Float64}(Iterators.flatten(local_results))
end

# Similarity of consensuses a and b. Upstream looks up (min, max) in a table keyed
# (max, min), so it always gets the default (deviation D6).
_similarity(table, a, b, default, upstream_quirks) =
    get(table, upstream_quirks ? (min(a, b), max(a, b)) : (max(a, b), min(a, b)), default)

# Upstream `calculate_match_lengths`: matched bases from the left and from the
# right of the alignment, each scan stopping after `chimera_allowable_errors`
# mismatches (not counting those within 15 bases of the query ends); a side is
# kept only if it is at least the minimum match length and longer than the other.
# 0-based half-open coordinates; returned in query orientation.
function _match_lengths(cigar, query, target, query_start, query_end, target_start, target_end, rc::Bool,
                        config::SavontConfig)
    allowed = config.chimera_allowable_errors
    left = 0; errors = 0; q = query_start; t = target_start
    for word in cigar
        errors > allowed && break
        n = Int(word >> 4); op = word & 0xf
        if op == 0
            for i in 0:(n - 1)
                if q + i < length(query) && t + i < length(target)
                    if query[q + i + 1] == target[t + i + 1]
                        left += 1
                    else
                        errors += 1
                        errors > allowed && q + i >= PCR_SLACK && break
                    end
                end
            end
            q += n; t += n
        elseif op == 1
            q += n
        elseif op == 2
            t += n
        else
            @warn "Unexpected CIGAR operation: $op"
        end
    end
    right = 0; errors = 0; q = query_end; t = target_end
    for word in Iterators.reverse(cigar)
        errors > allowed && break
        n = Int(word >> 4); op = word & 0xf
        if op == 0
            for i in 0:(n - 1)
                if query[q - i] == target[t - i]
                    right += 1
                else
                    errors += 1
                    errors > allowed && q - i + PCR_SLACK <= length(query) && break
                end
            end
            q -= n; t -= n
        elseif op == 1
            q -= n
        elseif op == 2
            t -= n
        else
            @warn "Unexpected CIGAR operation: $op"
        end
    end
    min_length = something(config.chimera_detect_length, max(config.min_read_length ÷ 10, 100))
    right_kept = right < min_length || left >= right ? nothing : right
    left_kept = left < min_length || right >= left ? nothing : left
    return rc ? (right_kept, left_kept) : (left_kept, right_kept)
end

"""
    detect_chimeras(consensuses, config) -> Vector{ChimeraInfo}

Upstream chimera detection (decompressed sequences). Each consensus is mapped
(`map-ont`) to every consensus more than 3 times deeper, and each mapping gives a
left or a right match length (`_match_lengths`). A consensus is chimeric when a
left match and a right match to two different parents cover at least
min(0.9 max(s, 0.7), 0.8) of it (and less than 1.5, or 1.8 if s < 0.99), where
the parents' similarity s is below 0.97, or below 0.995 if both are more than 10
times deeper; or when one match covers all but `chimera_detect_length` (default
`min_read_length` ÷ 10) bases of it and its estimated mismatches to that parent
exceed 20 / log2(depth ratio). Sets each consensus's `chimera_score` (the
smallest log(similarity) × match length). With `upstream_quirks`, similarities
are looked up as upstream, i.e. never found (deviation D6).
"""
function detect_chimeras(consensuses::AbstractVector{ConsensusSequence}, config::SavontConfig)
    isempty(consensuses) && return ChimeraInfo[]
    similarities = pairwise_similarities(consensuses)
    quirks = config.upstream_quirks
    similarity(a, b, default) = _similarity(similarities, a, b, default, quirks)
    min_match_length = something(config.chimera_detect_length, config.min_read_length ÷ 10)
    min_match_length < 10 &&
        @warn "Chimera detection match length is set to a very low value of $min_match_length < 10. This may lead to false positives."
    calls = [ChimeraInfo[] for _ in consensuses]
    scores = zeros(length(consensuses))
    Threads.@threads :dynamic for query in eachindex(consensuses)
        sequence = consensuses[query].decompressed
        depth = consensuses[query].depth
        len = length(sequence)
        lefts = Tuple{Int,Int}[]; rights = Tuple{Int,Int}[]      # (parent, match length)
        for parent in eachindex(consensuses)
            (parent == query || consensuses[parent].depth <= depth * 3) && continue
            reference = consensuses[parent].decompressed
            index, options = _seed_mapper(reference; bucket_bits=4)
            regions, _ = Minimap2.map_query(index, sequence, options)
            for r in regions
                (r.extra === nothing || isempty(r.extra.cigar)) && continue
                left, right = if r.rev
                    _match_lengths(r.extra.cigar, reverse_complement(sequence), reference, len - Int(r.qe), len - Int(r.qs),
                                   Int(r.rs), Int(r.re), true, config)
                else
                    _match_lengths(r.extra.cigar, sequence, reference, Int(r.qs), Int(r.qe), Int(r.rs), Int(r.re), false, config)
                end
                left === nothing || push!(lefts, (parent, left))
                right === nothing || push!(rights, (parent, right))
            end
        end
        score = 0.0
        for (parent, match) in Iterators.flatten((lefts, rights))
            s = similarity(parent, query, 1.0)
            s < 0.85 && match < 500 && continue
            score = min(score, log(s) * match)
        end
        scores[query] = score
        for (left_parent, left) in lefts
            for (right_parent, right) in rights
                left_parent == right_parent && continue
                s = similarity(left_parent, right_parent, 0.0)
                if s < 0.97 || (s < 0.995 && consensuses[left_parent].depth > depth * 10 && consensuses[right_parent].depth > depth * 10)
                    coverage = (left + right) / len
                    if coverage >= min(0.9 * max(s, 0.7), 0.8) && (coverage < 1.5 || (s < 0.99 && coverage < 1.8))
                        push!(calls[query], ChimeraInfo(query, left_parent, right_parent, left, right, len, coverage))
                        break
                    end
                end
            end
        end
        for (parent, match) in Iterators.flatten((lefts, rights))
            # usize subtraction upstream: a query shorter than the length never matches.
            len >= min_match_length && match >= len - min_match_length || continue
            mismatches = max(unsafe_trunc(Int, (1.0 - similarity(parent, query, 1.0)) * len), 0)
            ratio = consensuses[parent].depth / depth
            ratio < 3.0 && continue
            mismatches > 20.0 / log2(ratio) &&
                push!(calls[query], ChimeraInfo(query, parent, parent, match, 0, len, match / len))
        end
    end
    for (c, score) in zip(consensuses, scores)
        c.chimera_score = _as_i64(score)
    end
    return reduce(vcat, calls)
end

"""
    filter_chimeras(consensuses, chimeras) -> Vector{ConsensusSequence}

The consensuses that are not the query of any chimera call, in order.
"""
function filter_chimeras(consensuses::AbstractVector{ConsensusSequence}, chimeras::AbstractVector{ChimeraInfo})
    chimeric = Set(c.query for c in chimeras)
    return [c for (i, c) in enumerate(consensuses) if !(i in chimeric)]
end
