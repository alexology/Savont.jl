"""
    ReadRecord

One input read. `header` is the full header line, `identifier` its part up to
the first whitespace and `description` the rest. Reads whose last header field is
`rc` (Cutadapt's reverse-complement mark) are flagged by `reverse_complemented`
and, unless read with `cutadapt_rc=false`, stored reverse-complemented with
reversed qualities. `sample_index` is the
position of the input file and `ordinal` the 1-based position across all files.
"""
struct ReadRecord
    header::String
    identifier::String
    description::String
    sequence::Vector{UInt8}
    qualities::Union{Nothing,Vector{UInt8}}
    sample_index::Int
    ordinal::Int
    reverse_complemented::Bool
end

function _format_from_path(path::AbstractString)
    stem = endswith(lowercase(path), ".gz") ? splitext(path)[1] : path
    extension = lowercase(splitext(stem)[2])
    extension in (".fa", ".fasta", ".fna", ".ffn", ".faa", ".frn") && return :fasta
    extension in (".fq", ".fastq") && return :fastq
    throw(ArgumentError("cannot infer FASTA/FASTQ format from path: $path"))
end

# `header` is the full header line without the leading '>' or '@' (FASTX.description).
function _record(identifier_value, header_value, sequence_value, quality_value, sample_index, ordinal,
                 cutadapt_rc::Bool)
    identifier_string = String(identifier_value)
    header = String(header_value)
    description_string = String(strip(SubString(header, ncodeunits(identifier_string) + 1)))
    sequence_bytes = Vector{UInt8}(codeunits(String(sequence_value)))
    quality_bytes = quality_value === nothing ? nothing : Vector{UInt8}(codeunits(String(quality_value)))
    fields = split(header)
    is_rc = !isempty(fields) && last(fields) == "rc"
    if is_rc && cutadapt_rc
        sequence_bytes = reverse_complement(sequence_bytes)
        quality_bytes === nothing || reverse!(quality_bytes)
    end
    return ReadRecord(header, identifier_string, description_string, sequence_bytes, quality_bytes,
                      sample_index, ordinal, is_rc)
end

"""
    foreach_record(f, paths; cutadapt_rc=true)

Stream the reads of `paths` (FASTA/FASTQ, optionally gzip-compressed, format
inferred from the extension) and call `f(::ReadRecord)` on each, in file order,
without keeping them in memory. With `cutadapt_rc=false`, reads marked `rc` are
kept as written (upstream does so when building retained reads).
"""
function foreach_record(f, paths::AbstractVector{<:AbstractString}; cutadapt_rc::Bool=true)
    ordinal = 0
    for (sample_index, path) in enumerate(paths)
        format = _format_from_path(path)
        raw = open(path, "r")
        try
            stream = endswith(lowercase(path), ".gz") ? GzipDecompressorStream(raw) : raw
            reader = format === :fasta ? FASTX.FASTA.Reader(stream) : FASTX.FASTQ.Reader(stream)
            try
                for record in reader
                    ordinal += 1
                    qualities = record isa FASTX.FASTQ.Record ? FASTX.quality(String, record) : nothing
                    f(_record(FASTX.identifier(record), FASTX.description(record),
                              FASTX.sequence(String, record), qualities, sample_index, ordinal, cutadapt_rc))
                end
            finally
                close(reader)
            end
        finally
            close(raw)
        end
    end
    return nothing
end
foreach_record(f, path::AbstractString; kwargs...) = foreach_record(f, [path]; kwargs...)

"""
    read_records(paths) -> Vector{ReadRecord}

Read all records of `paths` into memory; see [`foreach_record`](@ref).
"""
function read_records(paths::AbstractVector{<:AbstractString})
    records = ReadRecord[]
    foreach_record(record -> push!(records, record), paths)
    return records
end
read_records(path::AbstractString) = read_records([path])
