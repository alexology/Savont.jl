# Reference databases and taxonomy tables shared by classify and sintax
# (upstream `taxonomy.rs` and `databases.rs`).

"""
    TaxonomyEntry

A lineage: the `tax_id` and the ranks of upstream's `TaxonomyEntry`. Also used for
the assignment of an ASV (upstream `TaxonomyAssignment`, same fields), where ranks
below the reached level hold `UNCLASSIFIED`.
"""
struct TaxonomyEntry
    tax_id::String
    species::String
    genus::String
    family::String
    order::String
    class::String
    phylum::String
    clade::String
    superkingdom::String
    subspecies::String
    species_subgroup::String
    species_group::String
end

"""
    TaxonomyDatabase

A reference database directory: its `keyword` (`emu-1`, `silva-138.2`,
`greengenes2-2024.09` or `unite-10.0`), the reference FASTA and the taxonomy
keyed as [`database_key`](@ref) extracts keys from reference names.
"""
struct TaxonomyDatabase
    keyword::String
    fasta_path::String
    taxonomy::Dict{String,TaxonomyEntry}
end

const DATABASE_MARKER = ".savont_db"

# keyword => (description, loader); the download definitions are in download.jl.
const DATABASE_KEYWORDS = ["emu-1", "silva-138.2", "greengenes2-2024.09", "unite-10.0"]

_unknown_keyword(keyword) = ArgumentError("Unknown database keyword '$keyword'. Available: $(join(DATABASE_KEYWORDS, ", "))")

"""
    database_key(db, header) -> Union{Nothing,String}

Taxonomy key of a reference sequence from its name or header line: the text
before the first `:` (EMU), the accession before the first `.` of the first word
(SILVA), or the whole trimmed header (GreenGenes2, UNITE).
"""
function database_key(db::TaxonomyDatabase, header::AbstractString)
    header = lstrip(==('>'), header)
    if db.keyword == "emu-1"
        return String(first(split(header, ':')))
    elseif db.keyword == "silva-138.2"
        words = split(header)
        isempty(words) && return nothing
        return String(first(split(first(words), '.')))
    end
    header = strip(header)
    return isempty(header) ? nothing : String(header)
end

"""
    load_database(directory; upstream_quirks=false) -> TaxonomyDatabase

Load a database directory. Its type comes from the `.savont_db` marker written by
[`download_database`](@ref) or, without one, from the directory name. A SILVA
taxonomy map compressed with gzip is read; upstream accepts the name but fails
on the content (deviation D10), and so does `upstream_quirks=true`.
"""
function load_database(directory::AbstractString; upstream_quirks::Bool=false)
    marker = joinpath(directory, DATABASE_MARKER)
    keyword = isfile(marker) ? String(strip(read(marker, String))) : basename(rstrip(normpath(directory), ['/', '\\']))
    keyword == "emu-1" && return _load_emu(directory)
    keyword == "silva-138.2" && return _load_silva(directory; upstream_quirks)
    keyword == "greengenes2-2024.09" && return _load_gg2(directory)
    keyword == "unite-10.0" && return _load_unite(directory)
    throw(_unknown_keyword(keyword))
end

_open_text(f, path) = open(path) do raw
    f(endswith(path, ".gz") ? GzipDecompressorStream(raw) : raw)
end

_require_file(path, what) = isfile(path) || throw(SystemError("$what not found: $path", 2))

function _first_file(predicate, directory, what)
    for name in readdir(directory)
        path = joinpath(directory, name)
        isfile(path) && predicate(name) && return path
    end
    throw(SystemError("No $what found in $directory", 2))
end

function _load_emu(directory)
    fasta = joinpath(directory, "species_taxid.fasta"); taxonomy_path = joinpath(directory, "taxonomy.tsv")
    _require_file(fasta, "FASTA file"); _require_file(taxonomy_path, "Taxonomy file")
    taxonomy = Dict{String,TaxonomyEntry}()
    open(taxonomy_path) do io
        for (number, line) in enumerate(eachline(io))
            number == 1 && continue
            fields = split(line, '\t')
            if length(fields) < 12
                @warn "Skipping malformed line $number: insufficient fields"
                continue
            end
            entry = TaxonomyEntry(String.(fields[1:12])...)
            taxonomy[entry.tax_id] = entry
        end
    end
    return TaxonomyDatabase("emu-1", fasta, taxonomy)
end

# SILVA's taxonomy map, plain or gzip-compressed (upstream accepts a `.txt.gz`
# name but reads it as plain text, deviation D10).
function _load_silva(directory; upstream_quirks::Bool=false)
    fasta = _first_file(n -> endswith(n, ".fasta") || endswith(n, ".fasta.gz") || endswith(n, ".fa.gz"), directory, "FASTA file")
    taxmap = _first_file(n -> startswith(n, "taxmap_") && (endswith(n, ".txt") || endswith(n, ".txt.gz")), directory, "taxmap file")
    upstream_quirks && endswith(taxmap, ".gz") && throw(ArgumentError("stream did not contain valid UTF-8: $taxmap"))
    taxonomy = Dict{String,TaxonomyEntry}()
    _open_text(taxmap) do io
        for (number, line) in enumerate(eachline(io))
            number == 1 && continue
            fields = split(line, '\t')
            if length(fields) < 6
                @warn "Skipping malformed Silva line $number: insufficient fields"
                continue
            end
            levels = split(fields[4], ';')
            level(i) = String(strip(i <= length(levels) ? levels[i] : "UNKNOWN"))
            taxonomy[String(fields[1])] = TaxonomyEntry(String(fields[6]), String(fields[5]), level(6), level(5), level(4),
                                                        level(3), level(2), "", level(1), "", "", "")
        end
    end
    return TaxonomyDatabase("silva-138.2", fasta, taxonomy)
end

# Ranks of a lineage with `x__` prefixes (the last occurrence of a rank wins);
# `kingdom` is the letter of the top rank (`d` for GreenGenes2, `k` for UNITE),
# returned under `'d'`.
function _prefixed_ranks(text, kingdom::Char; underscores=false)
    ranks = Dict{Char,String}()
    for level in split(text, ';')
        level = strip(level)
        for tag in (kingdom, 'p', 'c', 'o', 'f', 'g', 's')
            prefix = string(tag, "__")
            if startswith(level, prefix)
                value = String(chopprefix(level, prefix))
                ranks[tag == kingdom ? 'd' : tag] = underscores ? replace(value, '_' => ' ') : value
                break
            end
        end
    end
    return tag -> get(ranks, tag, "")
end

function _load_gg2(directory)
    fasta = _first_file(n -> endswith(n, ".fa.gz") || endswith(n, ".fasta.gz") || endswith(n, ".fa"), directory, ".fa.gz file")
    taxonomy = Dict{String,TaxonomyEntry}()
    fill(s) = isempty(s) ? "Greengenes_unannotated" : s
    _open_text(fasta) do io
        for line in eachline(io)
            startswith(line, '>') || continue
            key = String(strip(line[2:end]))
            isempty(key) && continue
            rank = _prefixed_ranks(key, 'd')
            genus = rank('g'); epithet = rank('s')
            species = !isempty(genus) && !isempty(epithet) ? "$genus $epithet" : epithet
            taxonomy[key] = TaxonomyEntry(key, fill(species), fill(genus), fill(rank('f')), fill(rank('o')), fill(rank('c')),
                                          fill(rank('p')), "", fill(rank('d')), "", "", "")
        end
    end
    return TaxonomyDatabase("greengenes2-2024.09", fasta, taxonomy)
end

function _load_unite(directory)
    fasta = joinpath(directory, "unite_sequences.fasta")
    isfile(fasta) || throw(SystemError("unite_sequences.fasta not found in $directory. Run download_database with \"unite-10.0\" to set up the database.", 2))
    taxonomy = Dict{String,TaxonomyEntry}()
    fill(s) = isempty(s) ? "UNITE_unannotated" : s
    open(joinpath(directory, "unite_taxonomy.tsv")) do io
        for (number, line) in enumerate(eachline(io))
            number == 1 && continue
            fields = split(line, '\t'; limit=2)
            sh_id = String(fields[1]); isempty(sh_id) && continue
            rank = _prefixed_ranks(length(fields) > 1 ? fields[2] : "", 'k'; underscores=true)
            taxonomy[sh_id] = TaxonomyEntry(sh_id, fill(rank('s')), fill(rank('g')), fill(rank('f')), fill(rank('o')),
                                            fill(rank('c')), fill(rank('p')), "", fill(rank('d')), "", "", "")
        end
    end
    return TaxonomyDatabase("unite-10.0", fasta, taxonomy)
end

"""
    assign_taxonomy(entry, identity, species_threshold, genus_threshold, asv_header, detailed_unclassified)

Upstream `TaxonomyAssignment::from_taxonomy_entry`: the ranks of `entry` down to
species (identity ≥ `species_threshold`), genus (≥ `genus_threshold`), family
(≥ 86.5), order (≥ 82), class (≥ 78.5) or phylum (≥ 75); lower ranks are
`UNCLASSIFIED` (`UNCLASSIFIED-(asv_header)` with `detailed_unclassified`), and
everything below 75% identity.
"""
function assign_taxonomy(entry::TaxonomyEntry, identity::Real, species_threshold::Real, genus_threshold::Real,
                         asv_header::AbstractString, detailed_unclassified::Bool)
    u = detailed_unclassified ? "UNCLASSIFIED-($asv_header)" : "UNCLASSIFIED"
    e = entry
    identity >= species_threshold && return e
    level = identity >= genus_threshold ? 6 : identity >= 86.5 ? 5 : identity >= 82.0 ? 4 :
            identity >= 78.5 ? 3 : identity >= 75.0 ? 2 : 0
    level == 0 && return TaxonomyEntry(e.tax_id, u, u, u, u, u, u, u, u, "", "", "")
    rank(l, value) = level >= l ? value : u
    return TaxonomyEntry(e.tax_id, u, rank(6, e.genus), rank(5, e.family), rank(4, e.order), rank(3, e.class),
                         e.phylum, e.clade, e.superkingdom, "", "", "")
end

"""
    AsvClassification

Classification of one ASV (upstream `AsvClassification`): `taxonomy` is
`nothing` for unclassified ASVs; `identity`, `nm` and `best_hit_tax_id` come
from the chosen alignment (classify only).
"""
struct AsvClassification
    asv_id::String
    asv_header::String
    hit_reference_id::String
    abundance::Float64
    best_hit_tax_id::Union{Nothing,String}
    identity::Union{Nothing,Float64}
    nm::Union{Nothing,Int}
    taxonomy::Union{Nothing,TaxonomyEntry}
end

# Header fields up to the first whitespace and the depth token after the last `_`.
_header_token(header::AbstractString) = (words = split(header); isempty(words) ? String(header) : String(first(words)))
_depth_string(header::AbstractString) = String(last(split(_header_token(header), '_')))
_depth_sum(token) = sum((something(tryparse(Int, part), 0) for part in split(token, '-')); init=0)
# Upstream `extract_depths_from_headers`: per-sample depths are summed, at least 1.
_header_depth(header::AbstractString) = max(_depth_sum(_depth_string(header)), 1)

"""
    read_fasta_records(path) -> Vector{Tuple{String,Vector{UInt8}}}

`(header, sequence)` of every record of a FASTA file, plain or gzip-compressed;
the header is the whole line after `>`, as upstream reads it with needletail.
"""
function read_fasta_records(path::AbstractString)
    records = Tuple{String,Vector{UInt8}}[]
    open(path) do raw
        gzipped = !eof(raw) && peek(raw) == 0x1f
        stream = gzipped ? GzipDecompressorStream(raw) : raw
        reader = FASTX.FASTA.Reader(stream)
        try
            for record in reader
                push!(records, (String(FASTX.description(record)), Vector{UInt8}(codeunits(FASTX.sequence(String, record)))))
            end
        finally
            close(reader)
        end
    end
    return records
end

# Rust's `{}` for f64: the shortest round-trip digits, never in exponent notation.
function _rust_float(x::Float64)
    isnan(x) && return "NaN"
    isinf(x) && return x > 0 ? "inf" : "-inf"
    text = repr(x)
    negative = startswith(text, '-'); negative && (text = text[2:end])
    mantissa, exponent = occursin('e', text) ? split(text, 'e') : (text, "0")
    integer, fraction = split(mantissa, '.')
    digits = integer * fraction; point = length(integer) + parse(Int, exponent)
    result = point <= 0 ? "0." * "0"^(-point) * digits :
             point >= length(digits) ? digits * "0"^(point - length(digits)) :
             digits[1:point] * "." * digits[(point + 1):end]
    if occursin('.', result)
        result = rstrip(rstrip(result, '0'), '.')
    end
    result = lstrip(result, '0')
    (isempty(result) || startswith(result, '.')) && (result = "0" * result)
    return (negative ? "-" : "") * result
end

# Rows aggregated by key in first-seen order, then stably sorted by decreasing
# abundance. Upstream aggregates in a HashMap, so rows of equal abundance come out
# in an arbitrary order; this order is deterministic.
function _aggregate_by_abundance(classifications, key)
    order = String[]; rows = Dict{String,Tuple{TaxonomyEntry,Float64}}()
    for c in classifications
        c.taxonomy === nothing && continue
        k = key(c.taxonomy)
        if haskey(rows, k)
            rows[k] = (rows[k][1], rows[k][2] + c.abundance)
        else
            push!(order, k); rows[k] = (c.taxonomy, c.abundance)
        end
    end
    return sort!([rows[k] for k in order]; by=last, rev=true, alg=MergeSort)
end

_species_key(t) = join((t.species, t.genus, t.family, t.order, t.class, t.phylum, t.clade, t.superkingdom), '|')
_genus_key(t) = join((t.genus, t.family, t.order, t.class, t.phylum, t.clade, t.superkingdom), '|')

"""
    write_species_abundance(path, classifications)

`species_abundance.tsv`: the relative abundance of every species-level lineage.
"""
function write_species_abundance(path::AbstractString, classifications::AbstractVector{AsvClassification})
    open(path, "w") do io
        println(io, "abundance\tspecies\tgenus\tfamily\torder\tclass\tphylum\tclade\tsuperkingdom")
        for (t, abundance) in _aggregate_by_abundance(classifications, _species_key)
            println(io, join((_rust_float(abundance), t.species, t.genus, t.family, t.order, t.class, t.phylum, t.clade,
                              t.superkingdom), '\t'))
        end
    end
    return path
end

"""
    write_genus_abundance(path, classifications)

`genus_abundance.tsv`: the relative abundance of every genus-level lineage.
"""
function write_genus_abundance(path::AbstractString, classifications::AbstractVector{AsvClassification})
    open(path, "w") do io
        println(io, "abundance\tgenus\tfamily\torder\tclass\tphylum\tclade\tsuperkingdom")
        for (t, abundance) in _aggregate_by_abundance(classifications, _genus_key)
            println(io, join((_rust_float(abundance), t.genus, t.family, t.order, t.class, t.phylum, t.clade,
                              t.superkingdom), '\t'))
        end
    end
    return path
end

# Upstream `write_species_abundance_pooled` and `write_genus_abundance_pooled`:
# per-sample relative abundances, lineages in key order and then by decreasing
# total. The genus key leaves out the superkingdom, as upstream.
function _write_pooled_abundance(path, classifications, per_sample, sample_names, species::Bool)
    n = length(sample_names)
    totals = [foldl(+, (get(s, k, 0) for s in per_sample); init=0) for k in 1:n]
    rows = Dict{String,Tuple{TaxonomyEntry,Vector{Float64}}}()
    for c in classifications
        t = c.taxonomy; t === nothing && continue
        key = species ? _species_key(t) : join((t.genus, t.family, t.order, t.class, t.phylum, t.clade), '|')
        index = something(tryparse(Int, replace(c.asv_id, r"^(ASV_)+" => "")), 0)
        entry = get!(() -> (t, zeros(n)), rows, key)
        for k in 1:n
            depth = index < length(per_sample) ? get(per_sample[index + 1], k, 0) : 0
            totals[k] > 0 && (entry[2][k] += depth / totals[k])
        end
    end
    ordered = sort!([rows[k] for k in sort!(collect(keys(rows)))]; by=r -> foldl(+, r[2]; init=0.0), rev=true, alg=MergeSort)
    open(path, "w") do io
        print(io, species ? "species\t" : "", "genus\tfamily\torder\tclass\tphylum\tclade\tsuperkingdom")
        foreach(name -> print(io, '\t', name), sample_names); println(io)
        for (t, abundances) in ordered
            print(io, species ? "$(t.species)\t" : "", join((t.genus, t.family, t.order, t.class, t.phylum, t.clade, t.superkingdom), '\t'))
            foreach(a -> @printf(io, "\t%.6f", a), abundances); println(io)
        end
    end
    return path
end

"""
    write_asv_mappings(path, classifications)

`asv_mappings.tsv` of classify: one row per classification with its alignment
identity, number of mismatches, taxonomy key, lineage and reference name.
"""
function write_asv_mappings(path::AbstractString, classifications::AbstractVector{AsvClassification})
    open(path, "w") do io
        println(io, "asv_header\tdepth\talignment_identity\tnumber_mismatches\ttax_id\tspecies\tgenus\tfamily\torder\tclass\tphylum\tclade\tsuperkingdom\treference")
        for c in classifications
            depth = _depth_string(c.asv_header)
            if c.taxonomy !== nothing
                c.identity === nothing && continue
                t = c.taxonomy
                println(io, join((c.asv_header, depth, @sprintf("%.2f", c.identity), something(c.nm, 0),
                                  something(c.best_hit_tax_id, "NA"), t.species, t.genus, t.family, t.order, t.class,
                                  t.phylum, t.clade, t.superkingdom, c.hit_reference_id), '\t'))
            else
                println(io, c.asv_header, '\t', depth, "\tNA\tNA\tNA", "\tUNCLASSIFIED"^9)
            end
        end
    end
    return path
end
