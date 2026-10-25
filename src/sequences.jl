@inline function _complement(base::UInt8)
    base == UInt8('A') || base == UInt8('a') ? UInt8('T') :
    base == UInt8('T') || base == UInt8('t') ? UInt8('A') :
    base == UInt8('C') || base == UInt8('c') ? UInt8('G') :
    base == UInt8('G') || base == UInt8('g') ? UInt8('C') : UInt8('N')
end

function reverse_complement(sequence::AbstractVector{UInt8})
    result = Vector{UInt8}(undef, length(sequence))
    @inbounds for i in eachindex(sequence)
        result[i] = _complement(sequence[end - i + 1])
    end
    return result
end
reverse_complement(sequence::AbstractString) = String(reverse_complement(Vector{UInt8}(codeunits(sequence))))

"""
    homopolymer_compress(sequence, do_hpc=true; upstream_quirks=false) -> (compressed, lengths)

Homopolymer compression with runs capped at 255. Without `do_hpc` the sequence is
returned unchanged with lengths 1; upstream instead appends a copy of the first
base (its final-run push also runs in that branch), which `upstream_quirks=true`
reproduces (deviation D4 in `docs/deviations.md`).
"""
function homopolymer_compress(sequence::AbstractVector{UInt8}, do_hpc::Bool=true; upstream_quirks::Bool=false)
    isempty(sequence) && return (UInt8[], UInt8[])
    if !do_hpc
        upstream_quirks || return (collect(sequence), fill(UInt8(1), length(sequence)))
        return (push!(collect(sequence), sequence[1]), fill(UInt8(1), length(sequence) + 1))
    end
    compressed = UInt8[]
    lengths = UInt8[]
    sizehint!(compressed, length(sequence)); sizehint!(lengths, length(sequence))
    current = sequence[1]; run = UInt8(1)
    @inbounds for i in 2:length(sequence)
        if sequence[i] == current && run < typemax(UInt8)
            run += UInt8(1)
        else
            push!(compressed, current); push!(lengths, run)
            current = sequence[i]; run = UInt8(1)
        end
    end
    push!(compressed, current); push!(lengths, run)
    sizehint!(compressed, length(compressed)); sizehint!(lengths, length(lengths))
    return compressed, lengths
end

"""
    homopolymer_decompress(sequence, lengths)

Expand a homopolymer-compressed sequence. As upstream, a length mismatch logs a
warning and returns a copy of `sequence` unchanged.
"""
function homopolymer_decompress(sequence::AbstractVector{UInt8}, lengths::AbstractVector{UInt8})
    if length(sequence) != length(lengths)
        @warn "HPC sequence and lengths mismatch: $(length(sequence)) vs $(length(lengths))"
        return collect(sequence)
    end
    result = UInt8[]
    sizehint!(result, sum(Int, lengths))
    for (base, run) in zip(sequence, lengths)
        append!(result, Iterators.repeated(base, Int(run)))
    end
    return result
end

"""
    homopolymer_compress_with_quality(sequence, qualities, do_hpc=true)

Homopolymer compression keeping the minimum quality of each run (runs capped at
255). As upstream, an empty sequence or a quality vector of different length
returns three empty vectors.
"""
function homopolymer_compress_with_quality(sequence::AbstractVector{UInt8}, qualities::AbstractVector{UInt8}, do_hpc::Bool=true)
    isempty(sequence) || length(sequence) != length(qualities) ? (return UInt8[], UInt8[], UInt8[]) : nothing
    do_hpc || return (collect(sequence), collect(qualities), fill(UInt8(1), length(sequence)))
    compressed = UInt8[]; compressed_quality = UInt8[]; lengths = UInt8[]
    current = sequence[1]; run = UInt8(1); minimum_quality = qualities[1]
    @inbounds for i in 2:length(sequence)
        if sequence[i] == current && run < typemax(UInt8)
            run += UInt8(1); minimum_quality = min(minimum_quality, qualities[i])
        else
            push!(compressed, current); push!(compressed_quality, minimum_quality); push!(lengths, run)
            current = sequence[i]; run = UInt8(1); minimum_quality = qualities[i]
        end
    end
    push!(compressed, current); push!(compressed_quality, minimum_quality); push!(lengths, run)
    return compressed, compressed_quality, lengths
end

function quality_pool(qualities::AbstractVector{UInt8})
    pooled = collect(qualities)
    length(qualities) < 5 && return pooled
    @inbounds for i in eachindex(qualities)
        # Exact translation of Rust's zero-based condition and half-open range.
        if i >= 4 && i <= length(qualities) - 2
            pooled[i] = minimum(@view qualities[i-2:i+1])
        end
    end
    return pooled
end

function estimate_sequence_identity(qualities::AbstractVector{UInt8})
    isempty(qualities) && return nothing
    total = 0.0
    @inbounds for quality in qualities
        total += 10.0 ^ (-((Int(quality) - 33) / 10.0))
    end
    return 100.0 - total / length(qualities) * 100.0
end
