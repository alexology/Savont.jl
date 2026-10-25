# Stage 4b: pileups of the reads on each consensus, quality error rates and
# posterior analysis of every position (upstream `alignment.rs`:
# generate_consensus_pileups, estimate_quality_error_rates,
# analyze_pileup_consensuses, collect_low_quality_consensuses).

const MAX_SEQS_CONSENSUS = 250
const MAX_INSERTION_LENGTH = 2
const DEFAULT_ERR_RATE = 0.02
const INDEL_QUALITY = UInt8(48)          # quality whose error rate scores deletions
const BAD_LENGTH_THRESHOLD = 100

const PILEUP_BASE = 0x00
const PILEUP_DELETION = 0x01
const PILEUP_INSERTION = 0x02

"""
    PileupEntry

One read at one consensus position: an aligned base (`kind == PILEUP_BASE`) with
its quality and homopolymer length, a deletion, or an insertion after the
position. An insertion keeps its first `inserted` bases (at most 2, as upstream):
the first in `base`, `quality`, `hp_length`, the second in the `*2` fields.
"""
struct PileupEntry
    kind::UInt8
    base::UInt8
    quality::UInt8
    hp_length::UInt8
    inserted::UInt8
    base2::UInt8
    quality2::UInt8
    hp_length2::UInt8
end
_pileup_base(base, quality, hp) = PileupEntry(PILEUP_BASE, base, quality, hp, 0, 0, 0, 0)
const _PILEUP_DELETION = PileupEntry(PILEUP_DELETION, 0, 0, 0, 0, 0, 0, 0)

"""
    Pileup

Reads aligned at one position of a consensus: `position` (1-based), the
consensus base, the median homopolymer length of the aligned bases (1 if none),
the entries in read order and, for low-confidence positions, the log posterior of
a base other than the consensus one.
"""
mutable struct Pileup
    position::Int
    ref_base::UInt8
    ref_hp_length::UInt8
    entries::Vector{PileupEntry}
    alt_posterior::Union{Nothing,Float64}
end

"""
    consensus_pileup(reads, consensus, config) -> Vector{Pileup}

Map the first 250 reads of the cluster (in cluster order) to the consensus with
minimap2 as upstream, and record, along the CIGAR of each read's first mapping,
its bases, deletions and insertions (first two inserted bases, attached to the
preceding position). Reads are homopolymer-compressed with `use_hpc`. Each
position's `ref_hp_length` is the median homopolymer length of its bases.
"""
function consensus_pileup(reads::AbstractVector{TwinRead}, consensus::ConsensusSequence, config::SavontConfig)
    cluster = consensus.cluster
    n_align = min(length(cluster), MAX_SEQS_CONSENSUS)
    reference = consensus.sequence
    len = length(reference)
    pileup = [Pileup(p, reference[p], consensus.hp_lengths[p], sizehint!(PileupEntry[], n_align), nothing) for p in 1:len]
    index, options = _seed_mapper(reference)
    for k in 1:n_align
        read = reads[cluster[k]]
        _, expanded = _expanded_qualities(read)
        sequence, quality, hp = homopolymer_compress_with_quality(read.sequence, expanded, config.use_hpc)
        regions, _ = Minimap2.map_query(index, sequence, options)
        isempty(regions) && continue
        best = regions[1]
        (best.extra === nothing || isempty(best.extra.cigar)) && continue
        if best.rev
            sequence = reverse_complement(sequence); quality = reverse(quality); hp = reverse(hp)
            query_start = length(sequence) - best.qe
        else
            query_start = Int(best.qs)
        end
        mapped_length = Int(best.qe - best.qs)
        ref_pos = Int(best.rs)          # 0-based, as upstream
        query_pos = 0                   # 0-based within the mapped part of the read
        for word in best.extra.cigar
            n = Int(word >> 4); op = word & 0xf
            if op == 0
                for j in 0:(n - 1)
                    if ref_pos + j < len && query_pos + j < mapped_length
                        q = query_start + query_pos + j + 1
                        push!(pileup[ref_pos + j + 1].entries, _pileup_base(sequence[q], quality[q], hp[q]))
                    end
                end
                ref_pos += n; query_pos += n
            elseif op == 1
                if ref_pos > 0 && ref_pos - 1 < len && query_pos + n <= mapped_length
                    q = query_start + query_pos + 1
                    entry = if n >= 2
                        PileupEntry(PILEUP_INSERTION, sequence[q], quality[q], hp[q], 2, sequence[q + 1], quality[q + 1], hp[q + 1])
                    else
                        PileupEntry(PILEUP_INSERTION, sequence[q], quality[q], hp[q], 1, 0, 0, 0)
                    end
                    push!(pileup[ref_pos].entries, entry)
                end
                query_pos += n
            elseif op == 2
                for j in 0:(n - 1)
                    ref_pos + j < len && push!(pileup[ref_pos + j + 1].entries, _PILEUP_DELETION)
                end
                ref_pos += n
            else
                @warn "Unexpected CIGAR operation in pileup: $op"
            end
        end
    end
    hp_lengths = UInt8[]
    for p in pileup
        empty!(hp_lengths)
        for e in p.entries
            e.kind == PILEUP_BASE && push!(hp_lengths, e.hp_length)
        end
        if isempty(hp_lengths)
            p.ref_hp_length = 0x01
        else
            sort!(hp_lengths); mid = length(hp_lengths) ÷ 2
            p.ref_hp_length = iseven(length(hp_lengths)) ?
                UInt8((UInt16(hp_lengths[mid]) + UInt16(hp_lengths[mid + 1])) ÷ 2) : hp_lengths[mid + 1]
        end
    end
    return pileup
end

"""
    top_cluster_indices(consensuses, fraction) -> Vector{Int}

Indices of the deepest `fraction` of consensuses (rounded half away from zero, as
Rust's `round`), by decreasing depth then index.
"""
function top_cluster_indices(consensuses::AbstractVector{ConsensusSequence}, fraction::Float64)
    order = sort(collect(eachindex(consensuses)); by=i -> (-consensuses[i].depth, i))
    return order[1:round(Int, fraction * length(consensuses), RoundNearestTiesAway)]
end

"""
    quality_error_rates(pileups) -> Dict{UInt8,Float64}

Empirical error rate per quality value, from positions with less than 5% errors
(mismatches, deletions and insertions); each quality starts from one error in
one base.
"""
function quality_error_rates(pileups)
    counts = Dict{UInt8,Tuple{Int,Int}}()
    foreach(pileup -> _quality_counts!(counts, pileup), pileups)
    return _error_rates(counts)
end

# Mismatches and bases per quality of one pileup, added to `counts`; integer
# counts, so pileups can be counted one at a time and in any order.
function _quality_counts!(counts::Dict{UInt8,Tuple{Int,Int}}, pileup)
    for p in pileup
        total = 0; errors = 0
        for e in p.entries
            total += 1
            (e.kind != PILEUP_BASE || e.base != p.ref_base) && (errors += 1)
        end
        total > 0 && errors / total < 0.05 || continue
        for e in p.entries
            e.kind == PILEUP_BASE || continue
            errors_q, total_q = get(counts, e.quality, (0, 0))
            counts[e.quality] = (errors_q + (e.base != p.ref_base), total_q + 1)
        end
    end
    return counts
end

# Each quality starts from one error in one base.
_error_rates(counts) = Dict(q => (e + 1) / (t + 1) for (q, (e, t)) in counts)

function _log_sum_exp(a::Float64, b::Float64)
    m = max(a, b)
    isinf(m) && m < 0 && return -Inf
    return m + log(exp(a - m) + exp(b - m))
end

"""
    analyze_pileup!(consensus, pileup, error_rates, config) -> pileup

Upstream's posterior analysis of one consensus. Ends whose depth is below
max(maximum depth / 3, 3/4 of `min_cluster_size`, 2) are trimmed. At each
remaining position, the log posterior of a base other than the consensus one is
computed from the pileup (bases weighted by the error rate of their quality,
deletions by that of quality 48, insertions by that of their first base);
positions above −min(`posterior_threshold_ln`, 3 `min_cluster_size`) are low
confidence. Low-confidence positions within 100 bases of the ends extend the
masked ends; the trimmed ends are masked with `N`; the other low-confidence
positions are recorded in `low_quality_positions` (and masked with
`mask_low_quality`). Returns the trimmed pileup.
"""
function analyze_pileup!(consensus::ConsensusSequence, pileup::Vector{Pileup}, error_rates::Dict{UInt8,Float64},
                         config::SavontConfig)
    isempty(pileup) && return pileup
    min_coverage_abs = max(config.min_cluster_size * 3 ÷ 4, 2)
    min_coverage = max(maximum(p -> length(p.entries), pileup) ÷ 3, min_coverage_abs)
    indel_error_rate = get(error_rates, INDEL_QUALITY, DEFAULT_ERR_RATE)
    rate(q) = get(error_rates, q, DEFAULT_ERR_RATE)
    # As upstream, without any position reaching min_coverage nothing is trimmed.
    first_kept = something(findfirst(p -> length(p.entries) >= min_coverage, pileup), 1)
    last_kept = something(findlast(p -> length(p.entries) >= min_coverage, pileup), length(pileup))
    pileup = pileup[first_kept:last_kept]
    threshold = min(config.posterior_threshold_ln, Float64(config.min_cluster_size * 3))
    for p in pileup
        log_ref = 0.0; log_not_ref = 0.0
        for e in p.entries
            if e.kind == PILEUP_BASE
                error_rate = rate(e.quality)
                accuracy = 1.0 - error_rate
                if e.base == p.ref_base
                    log_ref += log(accuracy); log_not_ref += log(error_rate)
                else
                    log_ref += log(error_rate); log_not_ref += log(accuracy)
                end
            elseif e.kind == PILEUP_DELETION
                log_ref += log(indel_error_rate); log_not_ref += log(1.0 - indel_error_rate)
            else
                # Only the base before the insertion counts (upstream takes 0 of the
                # inserted bases' own terms).
                error_rate = rate(e.quality)
                log_not_ref += log(1.0 - error_rate); log_ref += log(error_rate)
            end
        end
        alt_posterior = log_not_ref - _log_sum_exp(log_ref, log_not_ref)
        alt_posterior > -threshold && (p.alt_posterior = alt_posterior)
    end
    # Masking, in upstream's 0-based coordinates.
    low_confidence = [p.position - 1 for p in pileup if p.alt_posterior !== nothing]
    left_start = pileup[1].position - 1
    right_end = pileup[end].position
    start_polish = BAD_LENGTH_THRESHOLD + left_start
    # usize subtraction in upstream's release build: below 100 it wraps to a huge
    # value, so no position falls in the right end region.
    end_polish = right_end >= BAD_LENGTH_THRESHOLD ? right_end - BAD_LENGTH_THRESHOLD : typemax(Int)
    left_region = maximum((x for x in low_confidence if x < start_polish); init=-1)
    left_region == -1 && (left_region = left_start)
    right_region = minimum((x for x in low_confidence if x >= end_polish); init=typemax(Int))
    right_region == typemax(Int) && (right_region = right_end)
    sequence = consensus.sequence
    for pos in 0:(left_region - 1)
        sequence[pos + 1] = UInt8('N')
    end
    for pos in right_region:(length(sequence) - 1)
        sequence[pos + 1] = UInt8('N')
    end
    for x in low_confidence
        config.mask_low_quality && (sequence[x + 1] = UInt8('N'))
        left_region < x < right_region && push!(consensus.low_quality_positions, x + 1)
    end
    return pileup
end

"""
    is_low_quality(consensus, config)

Upstream's low-quality criterion: with `n` low-quality positions, n > 0 and
depth ÷ n² < `n_depth_cutoff`.
"""
function is_low_quality(consensus::ConsensusSequence, config::SavontConfig)
    n = length(consensus.low_quality_positions)
    return n > 0 && consensus.depth ÷ (n * n) < config.n_depth_cutoff
end

"""
    decompress!(consensus) -> consensus

Expand the homopolymer lengths into `decompressed`, without the leading and
trailing `N`s (an all-`N` sequence is kept whole, as upstream).
"""
function decompress!(consensus::ConsensusSequence)
    sequence = homopolymer_decompress(consensus.sequence, consensus.hp_lengths)
    first = something(findfirst(!=(UInt8('N')), sequence), 1)
    last = something(findlast(!=(UInt8('N')), sequence), length(sequence))
    consensus.decompressed = sequence[first:last]
    return consensus
end

"""
    polish_consensuses!(reads, consensuses, config) -> (kept, low_quality, error_rates)

Stage 4b. Quality error rates are estimated from the pileups of the deepest 10% of
consensuses; then every consensus gets its pileup's homopolymer lengths and the
posterior analysis (`analyze_pileup!`), and is decompressed. Consensuses meeting
`is_low_quality` are returned separately; both lists keep the input order.
Consensuses are processed in parallel on Julia's threads.
"""
function polish_consensuses!(reads::AbstractVector{TwinRead}, consensuses::AbstractVector{ConsensusSequence},
                             config::SavontConfig)
    # Error counts of the deepest 10% of the clusters, one pileup at a time per
    # task (upstream holds those pileups together).
    top = top_cluster_indices(consensuses, 0.1)
    partial = _tmap(top) do index
        _quality_counts!(Dict{UInt8,Tuple{Int,Int}}(), consensus_pileup(reads, consensuses[index], config))
    end
    counts = Dict{UInt8,Tuple{Int,Int}}()
    for part in partial, (q, (e, t)) in part
        errors_q, total_q = get(counts, q, (0, 0))
        counts[q] = (errors_q + e, total_q + t)
    end
    error_rates = _error_rates(counts)
    Threads.@threads :dynamic for consensus in consensuses
        pileup = consensus_pileup(reads, consensus, config)
        consensus.hp_lengths = UInt8[p.ref_hp_length for p in pileup]
        analyze_pileup!(consensus, pileup, error_rates, config)
    end
    foreach(decompress!, consensuses)
    low = [is_low_quality(c, config) for c in consensuses]
    return consensuses[.!low], consensuses[low], error_rates
end
