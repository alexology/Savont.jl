# `savont export` (upstream `merge.rs`): ASVs of several runs combined by
# sequence, with QIIME2-compatible tables.

# djb2 over upper-case bytes.
function _djb2(sequence::AbstractVector{UInt8})
    h = UInt64(5381)
    for b in sequence
        h = h * UInt64(33) + UInt64(UInt8('a') <= b <= UInt8('z') ? b - 0x20 : b)
    end
    return h
end

"""
    sequence_hash(sequence) -> String

Feature ID of an ASV in exported tables: the smaller djb2 hash (upper case) of
the sequence and of its reverse complement, as 16 hexadecimal digits.
"""
sequence_hash(sequence::AbstractVector{UInt8}) = string(min(_djb2(sequence), _djb2(reverse_complement(sequence))); base=16, pad=16)

# Feature table of a run directory: sample names and per-sample depths by OTU id.
function _export_feature_table(directory)
    path = joinpath(directory, "feature-table.tsv")
    isfile(path) || return nothing
    lines = readlines(path)
    start = findfirst(l -> startswith(l, "#OTU ID"), lines)
    start === nothing && return nothing
    names = String.(split(lines[start], '\t')[2:end])
    isempty(names) && return nothing
    depths = Dict{String,Vector{Int}}()
    for line in lines[(start + 1):end]
        (isempty(line) || startswith(line, '#')) && continue
        fields = split(line, '\t')
        depths[String(fields[1])] = [k + 1 <= length(fields) ? something(_parse_unsigned(fields[k + 1]), 0) : 0
                                     for k in eachindex(names)]
    end
    return names, depths
end

# `(header word, lineage)` of every row of asv_mappings.tsv (classify or sintax),
# the lineage joining the available ranks from superkingdom to species with `;`.
function _mapping_lineages(path)
    lines = readlines(path)
    isempty(lines) && throw(ArgumentError("empty file"))
    columns = split(lines[1], '\t')
    indices = [findfirst(==(name), columns) for name in ("superkingdom", "phylum", "class", "order", "family", "genus", "species")]
    pairs = Tuple{String,String}[]
    for line in lines[2:end]
        isempty(line) && continue
        fields = split(line, '\t')
        lineage = join((fields[i] for i in indices if i !== nothing && i <= length(fields)), ';')
        push!(pairs, (String(fields[1]), lineage))
    end
    return pairs
end

# Depths of an ASV missing from its run's feature table, from its FASTA header.
# Upstream puts the header total in the first sample column; with `per_sample`
# (a feature table with `width` columns exists, deviation D16) the
# dash-separated per-sample depths of a pooled header are used instead, and a
# header without one depth per column is an error.
function _header_counts(header, width, per_sample::Bool)
    token = _depth_string(header)
    (per_sample && width > 1) || return [_depth_sum(token)]
    parts = [_parse_unsigned(part) for part in split(token, '-')]
    length(parts) == width && all(!isnothing, parts) && return Int[p for p in parts]
    throw(ArgumentError("ASV header $(repr(_header_token(header))) is not in feature-table.tsv and does not give one depth per sample ($width)"))
end

_fuzzy_minimizers(sequence; upstream_quirks) = unique!(sort!(first(minimizer_seeds_positions(sequence, 28, 31; upstream_quirks))))

"""
    fuzzy_merge!(table, lineages; upstream_quirks=false) -> Int

Upstream `fuzzy_merge_table`: from the shortest ASV up, an ASV whose (w = 28,
k = 31) minimizers all occur in longer ASVs at most 10 bases longer (or of equal
length) is absorbed by the deepest of them: its counts are added and its lineage
passed on if the representative has none. `table` maps feature IDs to
`(sequence, counts)`, `lineages` feature IDs to lineages. Returns the number of
absorbed ASVs. Upstream picks an arbitrary representative among equally deep
candidates (hash-set order); here the smallest feature ID. With the corrected
minimizers (deviation D1) ASVs trimmed at the 5' end are merged too.
"""
function fuzzy_merge!(table::Dict{String,Tuple{Vector{UInt8},Vector{Int}}}, lineages::Dict{String,String};
                      upstream_quirks::Bool=false)
    ids = sort!(collect(keys(table)))
    minimizers = Dict(id => _fuzzy_minimizers(table[id][1]; upstream_quirks) for id in ids)
    inverted = Dict{UInt64,Set{String}}()
    for id in ids, kmer in minimizers[id]
        push!(get!(Set{String}, inverted, kmer), id)
    end
    absorbed = Set{String}()
    for id in sort(ids; by=id -> length(table[id][1]), alg=MergeSort)
        id in absorbed && continue
        kmers = minimizers[id]
        isempty(kmers) && continue
        len = length(table[id][1])
        candidates = copy(inverted[kmers[1]])
        for kmer in kmers
            intersect!(candidates, inverted[kmer])
            isempty(candidates) && break
        end
        delete!(candidates, id)
        filter!(c -> !(c in absorbed) && len <= length(table[c][1]) <= len + 10, candidates)
        isempty(candidates) && continue
        total(c) = foldl(+, table[c][2]; init=0)
        best = minimum(c for c in candidates if total(c) == maximum(total, candidates))
        table[best][2] .+= table[id][2]
        haskey(lineages, best) || !haskey(lineages, id) || (lineages[best] = lineages[id])
        for kmer in kmers
            delete!(inverted[kmer], id)
        end
        push!(absorbed, id)
    end
    for id in absorbed
        delete!(table, id); delete!(lineages, id)
    end
    if !isempty(absorbed)
        @info "Fuzzy merge: $(length(ids)) → $(length(table)) unique ASVs ($(length(absorbed)) near-identical sequences absorbed)"
    end
    return length(absorbed)
end

"""
    savont_export(input_dirs, output_dir; fuzzy=true, relabel=nothing, upstream_quirks=false)

Upstream `savont export`: combine the `final_asvs.fasta` of several runs (one
column per sample; a pooled run contributes all its columns) into ASVs identified
by [`sequence_hash`](@ref), so identical ASVs, also reverse complements, share a
row. Depths come from `feature-table.tsv`, or from the FASTA headers (for an ASV
missing from a pooled run's table, its per-sample header depths; deviation D16);
lineages
from `asv_mappings.tsv` of classify or sintax, the first one found for each ASV.
With `fuzzy`, near-identical ASVs are merged ([`fuzzy_merge!`](@ref)); `relabel`
replaces the sample names.

Writes `merged_feature_table.tsv`, `merged_rep_seqs.fasta`,
`merged_asv_taxonomy.tsv` and, when there is any ASV, `merged_taxon_counts.tsv`
(counts summed by lineage; an older one is removed otherwise). Returns `(table, lineages, sample_names)`.
"""
function savont_export(input_dirs::AbstractVector{<:AbstractString}, output_dir::AbstractString;
                       fuzzy::Bool=true, relabel::Union{Nothing,AbstractVector{<:AbstractString}}=nothing,
                       upstream_quirks::Bool=false)
    mkpath(output_dir)
    tables = [_export_feature_table(d) for d in input_dirs]
    sample_names = String[]; offsets = Int[]; widths = Int[]
    for (directory, ft) in zip(input_dirs, tables)
        push!(offsets, length(sample_names))
        names = ft === nothing ? [basename(rstrip(normpath(directory), ['/', '\\']))] : first(ft)
        push!(widths, length(names)); append!(sample_names, names)
    end
    n = length(sample_names)
    table = Dict{String,Tuple{Vector{UInt8},Vector{Int}}}()
    lineages = Dict{String,String}()
    for (d, directory) in enumerate(input_dirs)
        depths = tables[d] === nothing ? Dict{String,Vector{Int}}() : last(tables[d])
        fasta = joinpath(directory, ASV_FILE)
        records = try
            read_fasta_records(fasta)
        catch err
            @error "Could not read $fasta: $err"
            continue
        end
        token_hash = Dict{String,String}()
        unmatched = 0
        for (header, sequence) in records
            token = _header_token(header); id = sequence_hash(sequence)
            token_hash[token] = id
            counts = get(depths, token, nothing)
            if counts === nothing
                tables[d] === nothing || (unmatched += 1)
                counts = _header_counts(header, widths[d], tables[d] !== nothing && !upstream_quirks)
            end
            entry = get!(() -> (sequence, zeros(Int, n)), table, id)
            for (k, depth) in enumerate(counts)
                k <= widths[d] && (entry[2][offsets[d] + k] += depth)
            end
        end
        unmatched > 0 && @warn "$unmatched ASVs of $fasta are not in its feature-table.tsv; their depths are taken from the FASTA headers"
        mappings = joinpath(directory, "asv_mappings.tsv")
        if isfile(mappings)
            try
                for (token, lineage) in _mapping_lineages(mappings)
                    id = get(token_hash, token, nothing)
                    id === nothing || haskey(lineages, id) || (lineages[id] = lineage)
                end
            catch err
                @warn "Could not read $mappings: $err"
            end
        end
    end
    if relabel !== nothing
        length(relabel) == n || throw(ArgumentError("relabel: $(length(relabel)) label(s) provided for $n total sample column(s); counts must match"))
        sample_names = collect(String, relabel)
    end
    duplicates = sort!(unique!([s for (i, s) in enumerate(sample_names) if s in sample_names[1:(i - 1)]]))
    isempty(duplicates) || @warn "Duplicate sample names detected: [$(join(duplicates, ", "))]. Duplicate column names produce incorrect results in downstream tools; use `relabel` to assign unique names."
    fuzzy && fuzzy_merge!(table, lineages; upstream_quirks)

    ids = sort!(collect(keys(table)))
    open(joinpath(output_dir, "merged_feature_table.tsv"), "w") do io
        println(io, "#OTU ID", join(("\t" * s for s in sample_names)))
        foreach(id -> println(io, id, join(("\t$c" for c in table[id][2]))), ids)
    end
    open(joinpath(output_dir, "merged_rep_seqs.fasta"), "w") do io
        foreach(id -> write(io, '>', id, '\n', table[id][1], '\n'), ids)
    end
    open(joinpath(output_dir, "merged_asv_taxonomy.tsv"), "w") do io
        println(io, "Feature ID\tTaxon")
        foreach(id -> println(io, id, '\t', get(lineages, id, "Unclassified")), ids)
    end
    counts = Dict{String,Vector{Int}}()
    for id in ids
        counts[get(lineages, id, "Unclassified")] = get(() -> zeros(Int, n), counts, get(lineages, id, "Unclassified")) .+ table[id][2]
    end
    # Without ASVs there is no taxon table, as upstream; one left by an earlier
    # export is removed so that the files stay consistent (deviation D14).
    isempty(counts) && rm(joinpath(output_dir, "merged_taxon_counts.tsv"); force=true)
    if !isempty(counts)
        open(joinpath(output_dir, "merged_taxon_counts.tsv"), "w") do io
            println(io, "taxon", join(("\t" * s for s in sample_names)))
            foreach(l -> println(io, l, join(("\t$c" for c in counts[l]))), sort!(collect(keys(counts))))
        end
    end
    return (; table, lineages, sample_names)
end
