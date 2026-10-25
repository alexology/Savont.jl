# Stage 7: read assignment and EM refinement of ASV depths (upstream
# `alignment.rs`: refine_asv_depths_with_em, refine_asv_depths_with_minimap2).

"""
    ReadMapping

One row of upstream's `read_to_asv_mappings.tsv`: the `id` of an ASV the read was
assigned to, the SNPmer mismatches (`nothing` in the minimap2 path) and the edit
distance (NM) of the alignment. The first row of a read is its best assignment.
"""
struct ReadMapping
    asv_id::Int
    snpmer_mismatches::Union{Nothing,Int}
    nm::Int
end

"""
    asv_twin_reads(consensuses, info, config) -> Vector{TwinRead}

The decompressed consensuses as `TwinRead`s without qualities, with the
SNPmers of the stage-1 `info` and no blockmers (upstream `twin_reads_from_fasta`).
"""
function asv_twin_reads(consensuses::AbstractVector{ConsensusSequence}, info::KmerGlobalInfo, config::SavontConfig)
    snpmer_set = Set{UInt64}()
    for snpmer in info.snpmer_info
        union!(snpmer_set, snpmer_kmers(snpmer))
    end
    reads = TwinRead[]
    for (i, c) in enumerate(consensuses)
        read = twin_read(c.decompressed, nothing, "asv_$i", config, snpmer_set, Set{UInt64}())
        read === nothing && throw(ArgumentError("consensus $(c.id) is shorter than k"))
        push!(reads, read)
    end
    return reads
end

# Rust `f64::powi` (compiler-rt `__powidf2`): repeated squaring.
function _powi(a::Float64, n::Int)
    negative = n < 0; n = abs(n); r = 1.0
    while true
        isodd(n) && (r *= a)
        n ÷= 2
        n == 0 && break
        a *= a
    end
    return negative ? 1.0 / r : r
end

# EM over equivalence classes (sorted ASV index lists with read counts), as
# upstream: abundances start uniform, each round distributes every class's reads
# in proportion to the current abundances of its ASVs and divides by the
# assigned reads; it stops when no abundance moves by 0.01 / assigned reads, or
# after 10,000 rounds. Returns the depths, rounded so that they add up to the
# assigned reads (deviation D12); with `upstream_quirks`, each rounded on its own.
function _em_depths(classes::Vector{Pair{Vector{Int},Int}}, n::Int, assigned::Int; upstream_quirks::Bool=false)
    abundances = fill(1.0 / n, n)
    threshold = 0.01 / assigned
    for _ in 1:10_000
        updated = zeros(n)
        for (asvs, count) in classes
            denominator = 0.0
            for a in asvs
                denominator += abundances[a]
            end
            denominator > 0.0 || continue
            for a in asvs
                updated[a] += count * abundances[a] / denominator
            end
        end
        if sum(updated) > 0.0
            updated ./= assigned
        end
        change = 0.0
        for (old, new) in zip(abundances, updated)
            change = max(change, abs(old - new))
        end
        abundances = updated
        change < threshold && break
    end
    upstream_quirks && return [round(Int, a * assigned, RoundNearestTiesAway) for a in abundances]
    return _conserving_round(abundances .* assigned, assigned)
end

# Largest-remainder rounding: the integer parts, plus one for the largest
# fractional parts (earlier index on ties) until the total is reached.
function _conserving_round(values::AbstractVector{Float64}, total::Int)
    depths = [floor(Int, v) for v in values]
    remaining = total - sum(depths; init=0)
    if remaining > 0
        order = sortperm(values .- depths; rev=true, alg=MergeSort)
        for i in order[1:min(remaining, length(order))]
            depths[i] += 1
        end
    end
    return depths
end

# Best ASVs of one read by SNPmers, minimizers and alignment (upstream EM path):
# candidates sharing SNPmer contexts, kept if they share at least 0.95^k of the
# smaller minimizer set and have at most 0.005 SNPmer mismatches per shared
# minimizer and sampling step; among those with the fewest SNPmer mismatches,
# each ASV is mapped (lr:hq) to the read and the lowest NM wins.
function _snpmer_assignment(read::TwinRead, index, asv_minimizers, asv_sequences, k::Int, config::SavontConfig)
    mask = ~(UInt64(3) << (k - 1))
    statistics = Dict{Int,Tuple{Int,Int}}()
    for kmer in read.snpmer_kmers
        for (candidate, candidate_kmer) in get(index, kmer & mask, ())
            matches, mismatches = get(statistics, candidate, (0, 0))
            statistics[candidate] = kmer == candidate_kmer ? (matches + 1, mismatches) : (matches, mismatches + 1)
        end
    end
    minimizers = Set(read.minimizer_kmers)
    scores = Tuple{Int,Float64,Int}[]           # (asv, ratio, mismatches)
    for (asv, (_, mismatches)) in statistics
        shared = count(in(asv_minimizers[asv]), minimizers)
        shared == 0 && continue
        shared / min(length(minimizers), length(asv_minimizers[asv])) < _powi(0.95, k) && continue
        push!(scores, (asv, mismatches / shared / config.sampling_rate, mismatches))
    end
    best = [(asv, mismatches) for (asv, ratio, mismatches) in scores if ratio <= 0.005]
    isempty(best) && return ReadMapping[], Int[], typemax(Int32)
    fewest = minimum(last, best)
    filter!(x -> x[2] == fewest, best)
    index_options, map_options = Minimap2.preset("lr:hq")
    index_options.bucket_bits = 4
    read_index = Minimap2.MinimapIndex([read.sequence], index_options; names=["N/A"])
    map_options.mid_occ = 1000
    map_options.flag |= Minimap2.F_CIGAR
    alignments = Tuple{Int,Int,Int}[]          # (asv, nm, mismatches)
    for (asv, mismatches) in best
        regions, _ = Minimap2.map_query(read_index, asv_sequences[asv], map_options)
        isempty(regions) && continue
        r = regions[1]
        push!(alignments, (asv, Int(r.blen) - Int(r.mlen) + Int(r.extra.n_ambi), mismatches))
    end
    # Upstream orders ties by its hash-map iteration order, which decides the ASV
    # reported for reads tied between ASVs; here ties go to the lowest index.
    sort!(alignments; by=x -> (x[2], x[1]))
    best_nm = isempty(alignments) ? typemax(Int32) : alignments[1][2]
    asvs = sort!([asv for (asv, nm, _) in alignments if nm == best_nm])
    return [ReadMapping(asv, mismatches, nm) for (asv, nm, mismatches) in first(alignments, 5)], asvs, best_nm
end

# Best ASVs of one read by minimap2 (upstream low-polymorphism path): mappings
# with MAPQ > 0 to an index of all ASVs, ties at the lowest NM.
function _minimap2_assignment(read::TwinRead, index, options)
    regions, _ = Minimap2.map_query(index, read.sequence, options)
    valid = [r for r in regions if r.mapq > 0 && r.extra !== nothing]
    isempty(valid) && return ReadMapping[], Int[], typemax(Int32)
    nms = [Int(r.blen) - Int(r.mlen) + Int(r.extra.n_ambi) for r in valid]
    best_nm = minimum(nms)
    asvs = sort!(unique(Int(r.rid) + 1 for (r, nm) in zip(valid, nms) if nm == best_nm))
    return [ReadMapping(asv, nothing, best_nm) for asv in asvs], asvs, best_nm
end

# Best ASVs of every read, by the SNPmer path or, in low-polymorphism mode, by
# minimap2: per read, its mapping rows (ASV indices), best ASVs and best NM.
function _assign_reads(reads::AbstractVector{TwinRead}, consensuses::AbstractVector{ConsensusSequence},
                       info::KmerGlobalInfo, config::SavontConfig, low_polymorphism::Bool)
    n = length(consensuses)
    rows = Vector{Vector{ReadMapping}}(undef, length(reads))
    assignments = Vector{Vector{Int}}(undef, length(reads))
    nms = Vector{Int}(undef, length(reads))
    if low_polymorphism
        index_options, map_options = Minimap2.preset("lr:hq")
        index = Minimap2.MinimapIndex([c.decompressed for c in consensuses], index_options;
                                      names=["asv_$i" for i in 1:n])
        Minimap2.update_options!(map_options, index)
        map_options.flag |= Minimap2.F_CIGAR
        Threads.@threads :dynamic for r in eachindex(reads)
            rows[r], assignments[r], nms[r] = _minimap2_assignment(reads[r], index, map_options)
        end
    else
        asvs_as_reads = asv_twin_reads(consensuses, info, config)
        k = config.kmer_size; mask = ~(UInt64(3) << (k - 1))
        snpmer_index = Dict{UInt64,Vector{Tuple{Int,UInt64}}}()
        for (a, read) in enumerate(asvs_as_reads), (_, kmer) in snpmer_markers(read)
            push!(get!(Vector{Tuple{Int,UInt64}}, snpmer_index, kmer & mask), (a, kmer))
        end
        asv_minimizers = [Set(read.minimizer_kmers) for read in asvs_as_reads]
        asv_sequences = [read.sequence for read in asvs_as_reads]
        Threads.@threads :dynamic for r in eachindex(reads)
            rows[r], assignments[r], nms[r] = _snpmer_assignment(reads[r], snpmer_index, asv_minimizers, asv_sequences, k, config)
        end
    end
    return rows, assignments, nms
end

# Sorted equivalence classes (ASV index lists with read counts) of the assigned reads.
function _equivalence_classes(assignments)
    classes = Dict{Vector{Int},Int}()
    for asvs in assignments
        isempty(asvs) || (classes[asvs] = get(classes, asvs, 0) + 1)
    end
    return sort!(collect(classes); by=first)
end

"""
    refine_depths!(reads, consensuses, info, config; low_polymorphism=config.low_polymorphism)
        -> (consensuses, mappings)

Stage 7. Every retained read is assigned to its best ASVs, by SNPmers and
alignment (`_snpmer_assignment`), or in low-polymorphism mode by minimap2 against
all ASVs (`_minimap2_assignment`). Reads tied between several ASVs form
equivalence classes, and an EM estimate of the ASV abundances gives the new
depths; ASVs of depth 0 are removed. Also sets the counts of unambiguous,
ambiguous and NM ≤ 10 assignments. Returns the remaining consensuses (input
order) and, per read, its `ReadMapping` rows. Without any assigned read the
depths are kept. Reads are processed in parallel on Julia's threads.
"""
function refine_depths!(reads::AbstractVector{TwinRead}, consensuses::AbstractVector{ConsensusSequence},
                        info::KmerGlobalInfo, config::SavontConfig; low_polymorphism::Bool=config.low_polymorphism)
    isempty(consensuses) && return ConsensusSequence[], [ReadMapping[] for _ in reads]
    n = length(consensuses)
    rows, assignments, nms = _assign_reads(reads, consensuses, info, config, low_polymorphism)
    unambiguous = zeros(Int, n); ambiguous = zeros(Int, n); close = zeros(Int, n)
    for (asvs, nm) in zip(assignments, nms)
        isempty(asvs) && continue
        if length(asvs) == 1
            unambiguous[asvs[1]] += 1
        else
            ambiguous[asvs] .+= 1
        end
        nm <= 10 && (close[asvs] .+= 1)
    end
    for (i, c) in enumerate(consensuses)
        c.unambiguous_count = unambiguous[i]; c.ambiguous_count = ambiguous[i]; c.close_count = close[i]
    end
    # Rows refer to ASVs by id, as upstream's debug_id.
    mappings = [[ReadMapping(consensuses[m.asv_id].id, m.snpmer_mismatches, m.nm) for m in read_rows] for read_rows in rows]
    classes = _equivalence_classes(assignments)
    if isempty(classes)
        @warn "No reads mapped well to ASVs. Keeping original depths."
        return collect(consensuses), mappings
    end
    depths = _em_depths(classes, n, sum(last, classes); upstream_quirks=config.upstream_quirks)
    for (c, depth) in zip(consensuses, depths)
        c.depth = depth
    end
    return [c for c in consensuses if c.depth > 0], mappings
end

"""
    per_sample_depths(reads, consensuses, info, config, n_samples; low_polymorphism=config.low_polymorphism)
        -> Matrix{Int}

Stage 7b (`--pooled-samples`): the reads of each input file are assigned to the
final ASVs as in `refine_depths!` and an EM per sample gives the ASV × sample
depths (0 for samples without assigned reads).
"""
function per_sample_depths(reads::AbstractVector{TwinRead}, consensuses::AbstractVector{ConsensusSequence},
                           info::KmerGlobalInfo, config::SavontConfig, n_samples::Int;
                           low_polymorphism::Bool=config.low_polymorphism)
    n = length(consensuses)
    depths = zeros(Int, n, n_samples)
    (n == 0 || n_samples == 0) && return depths
    _, assignments, _ = _assign_reads(reads, consensuses, info, config, low_polymorphism)
    for sample in 1:n_samples
        classes = _equivalence_classes(a for (read, a) in zip(reads, assignments) if read.file_index == sample)
        isempty(classes) && continue
        depths[:, sample] = _em_depths(classes, n, sum(last, classes); upstream_quirks=config.upstream_quirks)
    end
    return depths
end
