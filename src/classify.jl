# `savont classify` (upstream `classify.rs`): ASVs aligned to a reference
# database with minimap2, multi-mapped ASVs resolved by EM over taxa, lineages cut
# by identity thresholds.

# One best-NM alignment of an ASV (upstream `AsvMapping`); `taxon` indexes the
# taxa in order of first appearance.
struct TaxonMapping
    asv::Int
    key::String
    taxon::Int
    reference::String
    identity::Float64
    nm::Int
    depth::Int
end

"""
    database_mapper(db) -> (index, options)

minimap2 index of the database FASTA, configured as upstream's
`Aligner::builder().map_ont().with_cigar().with_index(fasta)`: preset `map-ont`,
options updated from the index, base-level alignment. Reference names are the
header words before the first whitespace, as minimap2 reads them.
"""
function database_mapper(db::TaxonomyDatabase)
    records = read_fasta_records(db.fasta_path)
    isempty(records) && throw(ArgumentError("database FASTA is empty: $(db.fasta_path)"))
    index_options, map_options = Minimap2.preset("map-ont")
    index = Minimap2.MinimapIndex([s for (_, s) in records], index_options; names=[_header_token(h) for (h, _) in records])
    Minimap2.update_options!(map_options, index)
    map_options.flag |= Minimap2.F_CIGAR
    return index, map_options
end

_edit_distance(r::Minimap2.Region) = Int(r.blen) - Int(r.mlen) + Int(r.extra.n_ambi)

# Upstream `collect_best_mappings`: every alignment with the NM of the first
# (primary) one whose reference has a taxonomy entry, as
# (asv, key, identity, nm, reference name), ASVs in input order.
function _best_mappings(sequences, db, index, options)
    per_asv = _tmap(eachindex(sequences)) do asv
        found = Tuple{Int,String,Float64,Int,String}[]
        sequence = sequences[asv]
        isempty(sequence) && return found
        regions, _ = Minimap2.map_query(index, sequence, options)
        (isempty(regions) || regions[1].extra === nothing) && return found
        best = _edit_distance(regions[1])
        for r in regions
            r.extra === nothing && continue
            nm = _edit_distance(r)
            nm == best || continue
            identity = 100.0 * (1.0 - (nm / (r.qe - r.qs)))
            name = index.names[r.rid + 1]
            key = database_key(db, name)
            key !== nothing && haskey(db.taxonomy, key) && push!(found, (asv, key, identity, nm, name))
        end
        return found
    end
    return reduce(vcat, per_asv; init=Tuple{Int,String,Float64,Int,String}[])
end

"""
    taxon_em(mappings, n_taxa, total_reads, threshold) -> Vector{Float64}

Upstream `run_em_algorithm`: abundances of the taxa, starting uniform; each ASV
distributes its depth over its taxa in proportion to their current abundance,
the sums are divided by `total_reads`, until no abundance changes by
`threshold` or more (at most 1,000 iterations); abundances below `threshold` are
then set to 0. Upstream visits the ASVs in hash-map order; here they are
visited in input order.
"""
function taxon_em(mappings::AbstractVector{TaxonMapping}, n_taxa::Integer, total_reads::Integer, threshold::Real)
    abundances = fill(1.0 / n_taxa, n_taxa)
    groups = Dict{Int,Vector{Int}}()
    for (i, m) in enumerate(mappings)
        push!(get!(Vector{Int}, groups, m.asv), i)
    end
    order = sort!(collect(keys(groups)))
    iteration = 0
    while true
        iteration += 1
        updated = zeros(n_taxa)
        for asv in order
            members = groups[asv]
            denominator = foldl(+, (abundances[mappings[i].taxon] for i in members); init=0.0)
            denominator > 0 || continue
            for i in members
                m = mappings[i]
                updated[m.taxon] += m.depth * abundances[m.taxon] / denominator
            end
        end
        foldl(+, updated; init=0.0) > 0 && (updated ./= total_reads)
        change = foldl((x, y) -> max(x, y), (abs(a - b) for (a, b) in zip(abundances, updated)); init=0.0)
        abundances = updated
        (change < threshold || iteration >= 1000) && break
    end
    abundances[abundances .< threshold] .= 0.0
    return abundances
end

# Upstream `read_feature_table_for_classify`: sample names, per-ASV depths
# from feature-table.tsv, matched by the first header word (zeros for ASVs missing
# from the table), and the number of missing ASVs; `nothing` without a usable table.
function _read_feature_table(path, headers)
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
        depths[String(fields[1])] = [i + 1 <= length(fields) ? something(_parse_unsigned(fields[i + 1]), 0) : 0
                                     for i in eachindex(names)]
    end
    missing = count(h -> !haskey(depths, _header_token(h)), headers)
    return names, [get(() -> zeros(Int, length(names)), depths, _header_token(h)) for h in headers], missing
end

# Rust `str::parse::<usize>`: digits with an optional `+`, nothing else.
_parse_unsigned(text) = occursin(r"^\+?[0-9]+$", text) ? tryparse(Int, text) : nothing

"""
    savont_classify(input_dir, db; output_dir=input_dir, species_threshold=99.0,
                    genus_threshold=94.5, detailed_unclassified=false) -> NamedTuple

Upstream `savont classify`: align the ASVs of `input_dir/final_asvs.fasta` to the
database (a directory or a [`TaxonomyDatabase`](@ref)) with minimap2 (`map-ont`),
keep for each ASV the alignments with the NM of its primary alignment, resolve
multi-mapped ASVs by EM over taxa ([`taxon_em`](@ref)) and cut each lineage at
the rank its identity reaches ([`assign_taxonomy`](@ref)). Depths come from
`feature-table.tsv` when present, else from the FASTA headers.

Writes `species_abundance.tsv`, `genus_abundance.tsv` (one column per sample with
several samples) and `asv_mappings.tsv` (every best alignment, ASVs in input
order, alignments by decreasing EM abundance) to `output_dir`. Returns
`(classifications, mappings, abundances)`: one classification per ASV sorted by
decreasing abundance, the alignment rows, and the EM abundances of the taxa.

Ties between alignments of equal EM abundance go to the last one, as upstream.
Upstream caches the index as `<fasta>.mmi`; it is rebuilt here at every call.
ASVs missing from `feature-table.tsv` get depth 0 with a warning; if the total
depth is 0, an `ArgumentError` is thrown, while upstream (and
`upstream_quirks=true`) writes `NaN` abundances (deviation D15).
"""
function savont_classify(input_dir::AbstractString, db::TaxonomyDatabase; output_dir::AbstractString=input_dir,
                         species_threshold::Real=99.0, genus_threshold::Real=94.5, detailed_unclassified::Bool=false,
                         upstream_quirks::Bool=false)
    fasta = joinpath(input_dir, ASV_FILE)
    isfile(fasta) || throw(SystemError("Input FASTA not found: $fasta", 2))
    isdir(output_dir) || ispath(output_dir) || mkpath(output_dir)
    isdir(output_dir) || throw(ArgumentError("output path exists and is not a directory: $output_dir"))
    records = read_fasta_records(fasta)
    isempty(records) && @warn "No consensus sequences in $fasta"
    headers = first.(records)
    table = _read_feature_table(joinpath(input_dir, "feature-table.tsv"), headers)
    sample_names, per_sample, missing = table === nothing ? (["sample"], [[_header_depth(h)] for h in headers], 0) : table
    missing > 0 && @warn "$missing of $(length(headers)) ASVs of $fasta are not in feature-table.tsv; their depth is taken as 0"
    depths = [foldl(+, s; init=0) for s in per_sample]
    total_reads = foldl(+, depths; init=0)
    !isempty(records) && total_reads == 0 && !upstream_quirks &&
        throw(ArgumentError("total ASV depth is 0 in $(joinpath(input_dir, "feature-table.tsv")): abundances are undefined"))

    found = isempty(records) ? Tuple{Int,String,Float64,Int,String}[] :
            _best_mappings(last.(records), db, database_mapper(db)...)
    taxon_of = Dict{String,Int}(); taxa = String[]
    mappings = TaxonMapping[]
    for (asv, key, identity, nm, reference) in found
        taxon = get!(taxon_of, key) do
            push!(taxa, key); length(taxa)
        end
        push!(mappings, TaxonMapping(asv, key, taxon, reference, identity, nm, depths[asv]))
    end
    abundances = taxon_em(mappings, length(taxa), total_reads, 0.1 / total_reads)

    by_asv = [Int[] for _ in records]
    foreach(((i, m),) -> push!(by_asv[m.asv], i), enumerate(mappings))
    classifications = AsvClassification[]; rows = AsvClassification[]
    for (asv, (header, _)) in enumerate(records)
        asv_id = "ASV_$(asv - 1)"; asv_header = _header_token(header)
        abundance = depths[asv] / total_reads
        members = mappings[by_asv[asv]]
        if isempty(members)
            push!(classifications, AsvClassification(asv_id, asv_header, "", abundance, nothing, nothing, nothing, nothing))
            continue
        end
        classify(m) = AsvClassification(asv_id, asv_header, m.reference, abundance, m.key, m.identity, m.nm,
                                        assign_taxonomy(db.taxonomy[m.key], m.identity, species_threshold, genus_threshold,
                                                        asv_header, detailed_unclassified))
        append!(rows, classify.(sort(members; by=m -> abundances[m.taxon], rev=true, alg=MergeSort)))
        best = members[1]
        for m in members
            abundances[m.taxon] >= abundances[best.taxon] && (best = m)
        end
        push!(classifications, classify(best))
    end
    sort!(classifications; by=c -> c.abundance, rev=true, alg=MergeSort)

    if length(sample_names) > 1
        _write_pooled_abundance(joinpath(output_dir, "species_abundance.tsv"), classifications, per_sample, sample_names, true)
        _write_pooled_abundance(joinpath(output_dir, "genus_abundance.tsv"), classifications, per_sample, sample_names, false)
    else
        write_species_abundance(joinpath(output_dir, "species_abundance.tsv"), classifications)
        write_genus_abundance(joinpath(output_dir, "genus_abundance.tsv"), classifications)
    end
    write_asv_mappings(joinpath(output_dir, "asv_mappings.tsv"), rows)
    return (; classifications, mappings=rows, abundances)
end
savont_classify(input_dir::AbstractString, db::AbstractString; kwargs...) =
    savont_classify(input_dir, load_database(db); kwargs...)
