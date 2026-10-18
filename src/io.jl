struct ReadRecord
    identifier::String
    description::String
    sequence::Vector{UInt8}
    qualities::Union{Nothing,Vector{UInt8}}
    sample_index::Int
    ordinal::Int
    reverse_complemented::Bool
end

function _input_stream(path::AbstractString)
    raw = open(path, "r")
    endswith(lowercase(path), ".gz") ? GzipDecompressorStream(raw) : raw
end

function _format_from_path(path::AbstractString)
    stem = endswith(lowercase(path), ".gz") ? splitext(path)[1] : path
    extension = lowercase(splitext(stem)[2])
    extension in (".fa", ".fasta", ".fna", ".ffn", ".faa", ".frn") && return :fasta
    extension in (".fq", ".fastq") && return :fastq
    throw(ArgumentError("cannot infer FASTA/FASTQ format from path: $path"))
end

function _record(identifier_value, description_value, sequence_value, quality_value, sample_index, ordinal)
    identifier_string = String(identifier_value)
    description_string = String(description_value)
    sequence_bytes = Vector{UInt8}(codeunits(String(sequence_value)))
    quality_bytes = quality_value === nothing ? nothing : Vector{UInt8}(codeunits(String(quality_value)))
    header = isempty(description_string) ? identifier_string : string(identifier_string, ' ', description_string)
    is_rc = !isempty(header) && last(split(header)) == "rc"
    if is_rc
        sequence_bytes = reverse_complement(sequence_bytes)
        quality_bytes === nothing || reverse!(quality_bytes)
    end
    return ReadRecord(identifier_string, description_string, sequence_bytes, quality_bytes,
                      sample_index, ordinal, is_rc)
end

function read_records(paths::AbstractVector{<:AbstractString})
    records = ReadRecord[]
    ordinal = 0
    for (sample_index, path) in enumerate(paths)
        stream = _input_stream(path)
        reader = _format_from_path(path) === :fasta ? FASTX.FASTA.Reader(stream) : FASTX.FASTQ.Reader(stream)
        try
            for record in reader
                ordinal += 1
                qualities = record isa FASTX.FASTQ.Record ? FASTX.quality(String, record) : nothing
                push!(records, _record(FASTX.identifier(record), FASTX.description(record),
                                       FASTX.sequence(String, record), qualities,
                                       sample_index, ordinal))
            end
        finally
            close(reader)
        end
    end
    return records
end
read_records(path::AbstractString) = read_records([path])
