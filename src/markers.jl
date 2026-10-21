# Phase 3B: strand-separated k-mer counting and polymorphic marker inference
# (upstream `seq_parse::read_to_split_kmers`, `kmer_comp::get_snpmers_inplace_sort`
# and `kmer_comp::get_blockmers_inplace_sort`). Single-threaded reference path.

"""
    KmerCount

`(kmer, (count_a, count_b))` as produced by [`count_split_kmers`](@ref). For split
k-mers, `count_b` counts occurrences where the forward k-mer was canonical; for
blockmers, `count_b` counts forward-anchored occurrences.
"""
const KmerCount = Tuple{UInt64,NTuple{2,UInt32}}

"""
    SnpmerInfo

A biallelic SNPmer: the canonical k-mer with its middle base masked (`split_kmer`),
the two middle-base alleles ordered by decreasing abundance, their total counts and
the k-mer size.
"""
struct SnpmerInfo
    split_kmer::UInt64
    mid_bases::NTuple{2,UInt8}
    counts::NTuple{2,UInt32}
    k::UInt8
end
_order_key(x::SnpmerInfo) = (x.split_kmer, x.mid_bases, x.counts, x.k)
Base.isless(a::SnpmerInfo, b::SnpmerInfo) = isless(_order_key(a), _order_key(b))

"""
    snpmer_kmers(info::SnpmerInfo) -> NTuple{2,UInt64}

The two full (unmasked) k-mers of a SNPmer, as inserted into upstream's SNPmer set.
"""
snpmer_kmers(info::SnpmerInfo) = ntuple(i -> info.split_kmer | (UInt64(info.mid_bases[i]) << (info.k - 1)), 2)

Base.isless(a::Blockmer, b::Blockmer) = isless((a.kmer, a.is_forward), (b.kmer, b.is_forward))

"""
    BlockmerInfo

A biallelic blockmer: the anchor k-mer, the two full blockmers (anchor plus an
`l`-base suffix) ordered by decreasing abundance, their total counts, `k` and `l`.
"""
struct BlockmerInfo
    anchor_kmer::UInt64
    blockmers::NTuple{2,Blockmer}
    counts::NTuple{2,UInt32}
    k::UInt8
    l::UInt8
end
_order_key(x::BlockmerInfo) = (x.anchor_kmer, x.blockmers, x.counts, x.k, x.l)
Base.isless(a::BlockmerInfo, b::BlockmerInfo) = isless(_order_key(a), _order_key(b))

"""
    KmerGlobalInfo

Result of SNPmer inference: sorted `snpmer_info`, the set of high-frequency k-mers
(total count above `high_freq_threshold`) and the threshold itself.
"""
struct KmerGlobalInfo
    snpmer_info::Vector{SnpmerInfo}
    high_freq_kmers::Set{UInt64}
    high_freq_threshold::Float64
end

_total(count::KmerCount) = count[2][1] + count[2][2]

function _add_count!(index::Dict{UInt64,Int}, counts::Vector{KmerCount}, key::UInt64, second::Bool)
    slot = get!(index, key, length(counts) + 1)
    slot > length(counts) && push!(counts, (key, (UInt32(0), UInt32(0))))
    a, b = counts[slot][2]
    counts[slot] = (key, second ? (a, b + UInt32(1)) : (a + UInt32(1), b))
    return nothing
end

function _retain_counts(counts::Vector{KmerCount}, single_strand::Bool)
    keep = if single_strand
        filter(x -> x[2][1] > 2, counts)
    else
        filter(x -> x[2][1] > 0 && x[2][2] > 0 && _total(x) > 2, counts)
    end
    return sort!(keep; by=first)
end

"""
    count_split_kmers(paths, config::SavontConfig) -> (snpmer_counts, blockmer_counts)

Stream every read of `paths` once and count split k-mers (and, when
`config.use_blockmers`, blockmers) separately by orientation. Cutadapt `rc` reads
are reverse-complemented first. Counts are then filtered as upstream: with both
strands, a k-mer needs at least one occurrence per strand and more than two in
total; with `single_strand`, more than two occurrences in the first count.

Both results are `Vector{KmerCount}` sorted by k-mer. Bloom prefiltering
(`bloom_filter_size > 0`) is not implemented. Throws an `ErrorException` when
fewer than 0.1% of the split k-mers survive the filter, where upstream exits.
"""
function count_split_kmers(paths::AbstractVector{<:AbstractString}, config::SavontConfig)
    config.bloom_filter_size > 0 &&
        throw(ArgumentError("Bloom prefiltering (bloom_filter_size > 0) is not implemented"))
    k = config.kmer_size; l = config.blockmer_length; minimum_bq = config.minimum_base_quality
    mask = ~(UInt64(1) << 63)
    snpmer_index = Dict{UInt64,Int}(); snpmer_counts = KmerCount[]
    blockmer_index = Dict{UInt64,Int}(); blockmer_counts = KmerCount[]
    foreach_record(paths) do record
        for value in split_kmer_mid(record.sequence, record.qualities, k, minimum_bq)
            _add_count!(snpmer_index, snpmer_counts, value & mask, (value >> 63) == 1)
        end
        if config.use_blockmers
            for blockmer in blockmer_kmers(record.sequence, record.qualities, k, l, minimum_bq)
                _add_count!(blockmer_index, blockmer_counts, blockmer.kmer, blockmer.is_forward)
            end
        end
    end
    raw = length(snpmer_counts)
    snpmers = _retain_counts(snpmer_counts, config.single_strand)
    blockmers = _retain_counts(blockmer_counts, config.single_strand)
    if length(snpmers) < raw ÷ 1000
        error("less than 0.1% of split k-mers have counts on both strands and multiplicity > 2; " *
              "this may indicate a problem with the input data or very low coverage (consider single_strand=true)")
    end
    return snpmers, blockmers
end
count_split_kmers(path::AbstractString, config::SavontConfig) = count_split_kmers([path], config)

# Upstream sorts each group by decreasing total count, ties by k-mer, and keeps the
# top two alleles. Returns `nothing` when the pair is rejected.
function _biallelic_test(group::Vector{KmerCount}, binomial_short::Bool, single_strand::Bool)
    sort!(group; by=x -> (-Int64(_total(x)), x[1]))
    n = _total(group[1]); succ = _total(group[2])
    binomial_test(n, succ, 0.025) > 0.05 && return nothing
    binomial_short && binomial_test(n, succ, 0.050) > 0.05 && return nothing
    a = group[1][2][1]; b = group[2][2][1]; c = group[1][2][2]; d = group[2][2][2]
    table = (max(a, c), max(b, d), min(c, a), min(d, b))
    p_value = fishers_exact(table)
    odds = any(iszero, table) ? 0.0 : (Float64(table[1]) * Float64(table[4])) / (Float64(table[2]) * Float64(table[3]))
    !single_strand && odds == 0.0 && return nothing
    p_value > 0.005 || (odds < 1.5 && odds > 1.0 / 1.5) || return nothing
    return (n, succ)
end

function _each_group(f, sorted::Vector{KmerCount}, key, include)
    current = nothing
    group = KmerCount[]
    for entry in sorted
        include(entry) || continue
        group_key = key(entry[1])
        if current != group_key
            length(group) > 1 && f(group)
            group = KmerCount[]
            current = group_key
        end
        push!(group, entry)
    end
    length(group) > 1 && f(group)
    return nothing
end

"""
    find_snpmers(kmer_counts, config::SavontConfig) -> KmerGlobalInfo

Infer biallelic SNPmers from split k-mer counts (upstream
`get_snpmers_inplace_sort`). The high-frequency threshold is the count at the
upper 1/100,000 quantile, at least 100. K-mers sharing a masked context are
grouped; the two most abundant alleles are kept when the minor allele is too
frequent for sequencing error (binomial tail, p = 0.025) and its strand
distribution is not biased (two-sided Fisher test p > 0.005, or odds ratio within
(1/1.5, 1.5)). With `no_snpmers`, only the high-frequency set is computed.
"""
function find_snpmers(kmer_counts::AbstractVector{KmerCount}, config::SavontConfig)
    isempty(kmer_counts) && throw(ArgumentError("no k-mers found"))
    k = config.kmer_size
    totals = sort!([_total(x) for x in kmer_counts])
    n = length(totals)
    threshold = max(totals[n - n ÷ 100000], UInt32(100))
    high_freq = Set{UInt64}(x[1] for x in kmer_counts if _total(x) > threshold)
    config.no_snpmers && return KmerGlobalInfo(SnpmerInfo[], high_freq, Float64(threshold))
    sorted = sort(kmer_counts; by=x -> (split_kmer(x[1], k), x[1]))
    include = config.single_strand ? (_ -> true) : (x -> x[2][1] != 0 && x[2][2] != 0)
    snpmers = SnpmerInfo[]
    _each_group(sorted, x -> retrieve_masked_kmer(x, k), include) do group
        result = _biallelic_test(group, k < 5, config.single_strand)
        result === nothing && return
        mid_bases = (UInt8(split_kmer(group[1][1], k)[2]), UInt8(split_kmer(group[2][1], k)[2]))
        push!(snpmers, SnpmerInfo(retrieve_masked_kmer(group[1][1], k), mid_bases, result, UInt8(k)))
    end
    return KmerGlobalInfo(sort!(snpmers), high_freq, Float64(threshold))
end

"""
    find_blockmers(blockmer_counts, snpmer_counts, config::SavontConfig) -> Vector{BlockmerInfo}

Infer biallelic blockmers (upstream `get_blockmers_inplace_sort`). Each allele
needs more than two counts on each strand (first count only with
`single_strand`) and is skipped when its anchor k-mer, looked up among the
retained split k-mer counts, is more than ten times as frequent. Alleles sharing
an anchor are tested as in [`find_snpmers`](@ref), with the additional binomial
test at p = 0.05 when `blockmer_length < 5`. The result is sorted (upstream order
depends on thread scheduling).
"""
function find_blockmers(blockmer_counts::AbstractVector{KmerCount}, snpmer_counts::AbstractVector{KmerCount},
                        config::SavontConfig)
    k = config.kmer_size; l = config.blockmer_length; single_strand = config.single_strand
    anchor_counts = Dict{UInt64,UInt32}(x[1] => _total(x) for x in snpmer_counts)
    sorted = sort(blockmer_counts; by=first)
    function include(entry)
        a, b = entry[2]
        (a <= 2 || (!single_strand && b <= 2)) && return false
        return get(anchor_counts, entry[1] >> (2l), UInt32(0)) <= UInt32(10) * (a + b)
    end
    blockmers = BlockmerInfo[]
    _each_group(sorted, x -> x >> (2l), include) do group
        result = _biallelic_test(group, l < 5, single_strand)
        result === nothing && return
        pair = (Blockmer(group[1][1], true), Blockmer(group[2][1], true))
        push!(blockmers, BlockmerInfo(group[1][1] >> (2l), pair, result, UInt8(k), UInt8(l)))
    end
    return sort!(blockmers)
end

"""
    discover_markers(paths, config::SavontConfig) -> (KmerGlobalInfo, Vector{BlockmerInfo})

Count split k-mers and infer SNPmers and, when `config.use_blockmers`, blockmers
(upstream `get_kmers_and_snpmers`).
"""
function discover_markers(paths::AbstractVector{<:AbstractString}, config::SavontConfig)
    validate(config)
    snpmer_counts, blockmer_counts = count_split_kmers(paths, config)
    blockmers = config.use_blockmers ? find_blockmers(blockmer_counts, snpmer_counts, config) : BlockmerInfo[]
    return find_snpmers(snpmer_counts, config), blockmers
end
discover_markers(path::AbstractString, config::SavontConfig) = discover_markers([path], config)
