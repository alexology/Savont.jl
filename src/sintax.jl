# `savont sintax` (upstream `sintax.rs`): genus-level classification by k-mer
# bootstrap (SINTAX) against the database sequences.

const SINTAX_K = 12
const SINTAX_SUBSAMPLE = 32

# Canonical 12-mers (2-bit, `min` of the forward and reverse-complement codes) at
# every position; bases other than C, G, T and U count as A.
sintax_kmers(sequence::AbstractVector{UInt8}) = _sintax_kmers!(UInt32[], sequence)
function _sintax_kmers!(kmers::Vector{UInt32}, sequence::AbstractVector{UInt8})
    mask = (UInt32(1) << (2SINTAX_K)) - UInt32(1); shift = 2(SINTAX_K - 1)
    empty!(kmers); sizehint!(kmers, max(length(sequence) - SINTAX_K + 1, 0))
    forward = UInt32(0); reverse = UInt32(0)
    for (i, base) in enumerate(sequence)
        code = UInt32(encode_base(base))
        forward = ((forward << 2) | code) & mask
        reverse = (reverse >> 2) | ((UInt32(3) - code) << shift)
        i >= SINTAX_K && push!(kmers, min(forward, reverse))
    end
    return kmers
end

# Upstream's xorshift64 generator, seeded with max(seed, 1).
mutable struct _Xorshift
    state::UInt64
end
_Xorshift(seed::Integer) = _Xorshift(max(UInt64(seed), UInt64(1)))
function _next!(rng::_Xorshift)
    x = rng.state
    x ⊻= x << 13; x ⊻= x >> 7; x ⊻= x << 17
    rng.state = x
    return x
end

"""
    SintaxHit

SINTAX result of one ASV: header word, depth (from the header), relative
abundance, and for each rank the most frequent name among the bootstrap
iterations with its frequency.
"""
struct SintaxHit
    asv_header::String
    depth::Int
    abundance::Float64
    names::NTuple{7,String}       # species, genus, family, order, class, phylum, superkingdom
    bootstraps::NTuple{7,Float64}
end

_ranks(e::TaxonomyEntry) = (e.species, e.genus, e.family, e.order, e.class, e.phylum, e.superkingdom)

# Most frequent name; upstream takes an arbitrary one among equally frequent
# names (hash-map order), here the first one seen.
function _top_vote(names)
    counts = Dict{String,Int}(); order = String[]
    for name in names
        haskey(counts, name) || push!(order, name)
        counts[name] = get(counts, name, 0) + 1
    end
    best = order[1]
    for name in order
        counts[name] > counts[best] && (best = name)
    end
    return best, counts[best]
end

"""
    sintax_classify(sequences, depths, db; n_iter=100) -> Vector{Union{Nothing,Tuple}}

SINTAX bootstrap of upstream `sintax.rs`: for each ASV and iteration, 32 of its
canonical 12-mers drawn with replacement by a xorshift generator seeded with
`asv · n_iter + iteration + 1`; the database reference sharing most of them
(counted with multiplicity, each reference k-mer once) gives that iteration's
lineage, the earliest reference on ties (upstream processes references in
parallel, so its ties depend on thread timing; with one thread it is the same).
Returns per ASV `nothing` (no iteration classified) or the per-rank winning names
and their frequencies among the `n_iter` iterations.
"""
function sintax_classify(sequences::AbstractVector{<:AbstractVector{UInt8}}, db::TaxonomyDatabase; n_iter::Integer=100)
    n_pairs = length(sequences) * n_iter
    # The drawn k-mers as (k-mer, 0-based pair index asv · n_iter + iteration),
    # then indexed by k-mer: the pairs of k-mer x are pairs[starts[x + 1]:(starts[x + 2] - 1)].
    drawn = Tuple{UInt32,Int32}[]
    for (asv, sequence) in enumerate(sequences)
        kmers = sintax_kmers(sequence)
        isempty(kmers) && continue
        for iteration in 0:(n_iter - 1)
            pair = (asv - 1) * n_iter + iteration
            rng = _Xorshift(pair + 1)
            for _ in 1:SINTAX_SUBSAMPLE
                push!(drawn, (kmers[_next!(rng) % UInt64(length(kmers)) + 1], Int32(pair)))
            end
        end
    end
    length(drawn) < typemax(Int32) || throw(ArgumentError("too many ASVs × iterations for SINTAX"))
    starts = zeros(Int32, 4^SINTAX_K + 1)
    foreach(((kmer, _),) -> starts[kmer + 2] += Int32(1), drawn)
    starts[1] = 1; cumsum!(starts, starts)
    pairs = Vector{Int32}(undef, length(drawn)); filled = copy(starts)
    for (kmer, pair) in drawn
        pairs[filled[kmer + 1]] = pair; filled[kmer + 1] += Int32(1)
    end
    references = read_fasta_records(db.fasta_path)
    entries = [(key = database_key(db, h); key === nothing ? nothing : get(db.taxonomy, key, nothing)) for (h, _) in references]
    chunks = collect(Iterators.partition(eachindex(references), max(cld(length(references), 4Threads.nthreads()), 1)))
    partial = _tmap(chunks) do chunk
        chunk_scores = zeros(UInt16, n_pairs); chunk_best = zeros(Int32, n_pairs)
        counts = zeros(UInt16, n_pairs); touched = Int32[]; kmers = UInt32[]
        for r in chunk
            entries[r] === nothing && continue
            # Each distinct reference k-mer once.
            sort!(_sintax_kmers!(kmers, references[r][2]))
            for (i, kmer) in enumerate(kmers)
                i > 1 && kmers[i - 1] == kmer && continue
                for j in starts[kmer + 1]:(starts[kmer + 2] - 1)
                    pair = pairs[j]
                    counts[pair + 1] == 0 && push!(touched, pair)
                    counts[pair + 1] += UInt16(1)
                end
            end
            for pair in touched
                if counts[pair + 1] > chunk_scores[pair + 1]
                    chunk_scores[pair + 1] = counts[pair + 1]; chunk_best[pair + 1] = r
                end
                counts[pair + 1] = 0
            end
            empty!(touched)
        end
        return chunk_scores, chunk_best
    end
    scores = zeros(UInt16, n_pairs); best = zeros(Int32, n_pairs)
    for (s, b) in partial, pair in 1:n_pairs
        if s[pair] > scores[pair]
            scores[pair] = s[pair]; best[pair] = b[pair]
        end
    end
    return map(eachindex(sequences)) do asv
        lineages = [_ranks(entries[best[p]]) for p in ((asv - 1) * n_iter + 1):(asv * n_iter) if scores[p] > 0]
        isempty(lineages) && return nothing
        votes = [_top_vote(getindex.(lineages, rank)) for rank in 1:7]
        return (Tuple(first.(votes)), Tuple(last.(votes) ./ n_iter))
    end
end

"""
    savont_sintax(input_dir, db; output_dir=input_dir, min_bootstrap=0.8, n_iter=100,
                  detailed_unclassified=false) -> Vector{SintaxHit}

Upstream `savont sintax`: genus-level SINTAX classification
([`sintax_classify`](@ref)) of `input_dir/final_asvs.fasta` against the database
(a directory or a [`TaxonomyDatabase`](@ref)). Writes `genus_abundance.tsv`
(relative abundances by depths in the FASTA headers; ranks with bootstrap below
`min_bootstrap` are `UNCLASSIFIED`) and `asv_mappings.tsv` (bootstrap values and
names of every classified ASV, by decreasing abundance) to `output_dir`. Returns
the classified ASVs in that order. Species are never assigned.
"""
function savont_sintax(input_dir::AbstractString, db::TaxonomyDatabase; output_dir::AbstractString=input_dir,
                       min_bootstrap::Real=0.8, n_iter::Integer=100, detailed_unclassified::Bool=false)
    fasta = joinpath(input_dir, ASV_FILE)
    isfile(fasta) || throw(SystemError("Input FASTA not found: $fasta", 2))
    records = read_fasta_records(fasta)
    mkpath(output_dir)
    # Without sequences the tables are written with their headers only, so that no
    # earlier result is left behind (upstream stops with an error; deviation D14).
    isempty(records) && @warn "No sequences in $fasta"
    depths = [_header_depth(h) for (h, _) in records]
    total_reads = foldl(+, depths; init=0)
    results = isempty(records) ? [] : sintax_classify(last.(records), db; n_iter)
    hits = SintaxHit[]
    for (asv, result) in enumerate(results)
        result === nothing && continue
        abundance = total_reads > 0 ? depths[asv] / total_reads : 0.0
        push!(hits, SintaxHit(_header_token(lstrip(==('>'), records[asv][1])), depths[asv], abundance, result...))
    end
    sort!(hits; by=h -> h.abundance, rev=true, alg=MergeSort)

    classifications = map(hits) do h
        u = detailed_unclassified ? "UNCLASSIFIED-($(h.asv_header))" : "UNCLASSIFIED"
        call(rank) = h.bootstraps[rank] >= min_bootstrap ? h.names[rank] : u
        taxonomy = TaxonomyEntry("", u, call(2), call(3), call(4), call(5), call(6), "", call(7), "", "", "")
        AsvClassification(h.asv_header, h.asv_header, "", h.abundance, nothing, nothing, nothing, taxonomy)
    end
    write_genus_abundance(joinpath(output_dir, "genus_abundance.tsv"), classifications)
    open(joinpath(output_dir, "asv_mappings.tsv"), "w") do io
        println(io, "asv_header\tdepth\tspecies_bootstrap\tgenus_bootstrap\tfamily_bootstrap\torder_bootstrap\tclass_bootstrap\tphylum_bootstrap\tsuperkingdom_bootstrap\tspecies\tgenus\tfamily\torder\tclass\tphylum\tsuperkingdom")
        for h in hits
            call(rank) = h.bootstraps[rank] >= min_bootstrap ? h.names[rank] : "UNCLASSIFIED"
            println(io, h.asv_header, '\t', h.depth, join((@sprintf("\t%.3f", b) for b in h.bootstraps)),
                    "\tUNCLASSIFIED", join(("\t" * call(rank) for rank in 2:7)))
        end
    end
    return hits
end
savont_sintax(input_dir::AbstractString, db::AbstractString; kwargs...) =
    savont_sintax(input_dir, load_database(db); kwargs...)
