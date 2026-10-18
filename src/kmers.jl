struct Blockmer
    kmer::UInt64
    is_forward::Bool
end

@inline function encode_base(base::UInt8)::UInt64
    base == 0x01 && return 0x01
    base == 0x02 && return 0x02
    base == 0x03 && return 0x03
    base == UInt8('C') || base == UInt8('c') ? UInt64(1) :
    base == UInt8('G') || base == UInt8('g') ? UInt64(2) :
    base == UInt8('T') || base == UInt8('t') || base == UInt8('U') || base == UInt8('u') ? UInt64(3) : UInt64(0)
end

function decode_kmer(kmer::UInt64, k::Integer)
    result = Vector{UInt8}(undef, k)
    alphabet = UInt8['A', 'C', 'G', 'T']
    @inbounds for i in 1:k
        shift = 2 * (k - i)
        result[i] = alphabet[Int((kmer >> shift) & UInt64(3)) + 1]
    end
    return String(result)
end

@inline function mm_hash64(kmer::UInt64)::UInt64
    key = ~kmer + (kmer << 21)
    key = key ⊻ (key >> 24)
    key = (key + (key << 3)) + (key << 8)
    key = key ⊻ (key >> 14)
    key = (key + (key << 2)) + (key << 4)
    key = key ⊻ (key >> 28)
    return key + (key << 31)
end

@inline function rev_hash64(hashed_key::UInt64)::UInt64
    key = hashed_key
    tmp = key - (key << 31); key = key - (tmp << 31)
    tmp = key ⊻ (key >> 28); key = key ⊻ (tmp >> 28)
    key *= UInt64(14933078535860113213)
    tmp = key ⊻ (key >> 14); tmp = key ⊻ (tmp >> 14); tmp = key ⊻ (tmp >> 14); key = key ⊻ (tmp >> 14)
    key *= UInt64(15244667743933553977)
    tmp = key ⊻ (key >> 24); key = key ⊻ (tmp >> 24)
    tmp = ~key; tmp = ~(key - (tmp << 21)); tmp = ~(key - (tmp << 21)); return ~(key - (tmp << 21))
end

@inline function _parameters(k::Int)
    0 < k <= 31 || throw(ArgumentError("k must be in 1:31"))
    mask = typemax(UInt64) >> (64 - 2k)
    return mask, ~(UInt64(3) << (2k - 2)), 2(k - 1)
end

function fmh_seeds(sequence::AbstractVector{UInt8}, c::Integer, k::Integer)
    c > 0 || throw(ArgumentError("c must be positive")); length(sequence) < k && return UInt64[]
    mask, reverse_mask, reverse_shift = _parameters(k)
    forward = UInt64(0); reverse = UInt64(0); result = UInt64[]
    threshold = typemax(UInt64) ÷ UInt64(c)
    @inbounds for i in eachindex(sequence)
        nucleotide = encode_base(sequence[i]); complement = UInt64(3) - nucleotide
        forward = ((forward << 2) | nucleotide) & mask
        reverse = ((reverse >> 2) & reverse_mask) | (complement << reverse_shift)
        if i >= k
            hashed = mm_hash64(min(forward, reverse))
            hashed < threshold && push!(result, hashed)
        end
    end
    return result
end
fmh_seeds(sequence::AbstractString, c::Integer, k::Integer) = fmh_seeds(Vector{UInt8}(codeunits(sequence)), c, k)

function fmh_seeds_positions(sequence::AbstractVector{UInt8}, c::Integer, k::Integer)
    c > 0 || throw(ArgumentError("c must be positive")); length(sequence) < k && return (UInt64[], UInt64[])
    mask, reverse_mask, reverse_shift = _parameters(k)
    forward = UInt64(0); reverse = UInt64(0); kmers = UInt64[]; positions = UInt64[]
    threshold = typemax(UInt64) ÷ UInt64(c)
    @inbounds for i in eachindex(sequence)
        nucleotide = encode_base(sequence[i]); complement = UInt64(3) - nucleotide
        forward = ((forward << 2) | nucleotide) & mask
        reverse = ((reverse >> 2) & reverse_mask) | (complement << reverse_shift)
        if i >= k
            canonical = min(forward, reverse)
            if mm_hash64(canonical) < threshold
                push!(kmers, canonical); push!(positions, UInt64(i - k))
            end
        end
    end
    return kmers, positions
end

@inline _last_minimum_index(values) = findlast(==(minimum(values)), values)
function minimizer_seeds_positions(sequence::AbstractVector{UInt8}, w::Integer, k::Integer)
    w > 0 || throw(ArgumentError("w must be positive")); length(sequence) < k + w - 1 && return (UInt64[], UInt64[])
    mask, reverse_mask, reverse_shift = _parameters(k)
    forward = UInt64(0); reverse = UInt64(0); canonical = UInt64(0)
    window = fill(typemax(UInt64), w)
    @inbounds for i in 1:(k + w - 1)
        nucleotide = encode_base(sequence[i]); complement = UInt64(3) - nucleotide
        forward = (forward << 2) | nucleotide
        reverse = (reverse >> 2) | (complement << reverse_shift)
        if i >= k
            canonical = min(forward, reverse)
            window[i + 1 - k] = mm_hash64(canonical)
        end
    end
    minimum_slot = _last_minimum_index(window); minimum_value = window[minimum_slot]
    kmers = UInt64[canonical] # Preserve upstream's first-window behavior.
    positions = UInt64[minimum_slot - 1]
    @inbounds for i in (k + w):length(sequence)
        nucleotide = encode_base(sequence[i]); complement = UInt64(3) - nucleotide
        forward = ((forward << 2) | nucleotide) & mask
        reverse = ((reverse >> 2) & reverse_mask) | (complement << reverse_shift)
        hashed = mm_hash64(min(forward, reverse)); global_position = i - k
        window[mod(global_position, w) + 1] = hashed
        if hashed < minimum_value
            minimum_value = hashed; minimum_slot = mod(global_position, w) + 1
            push!(kmers, hashed); push!(positions, UInt64(global_position))
        elseif minimum_slot == mod(global_position, w) + 1
            minimum_slot = _last_minimum_index(window); minimum_value = window[minimum_slot]
            offset = mod(mod(global_position, w) - (minimum_slot - 1), w)
            push!(kmers, minimum_value); push!(positions, UInt64(global_position - offset))
        end
    end
    return kmers, positions
end

@inline retrieve_masked_kmer(kmer::UInt64, k::Integer) = kmer & ~(UInt64(3) << (k - 1))
@inline split_kmer(kmer::UInt64, k::Integer) = (retrieve_masked_kmer(kmer, k), (kmer >> (k - 1)) & UInt64(3))

function split_kmer_mid(sequence::AbstractVector{UInt8}, qualities::Union{Nothing,AbstractVector{UInt8}}, k::Integer, minimum_bq::Integer)
    length(sequence) < k && return UInt64[]
    isodd(k) && k <= 31 || throw(ArgumentError("k must be odd and <= 31"))
    qualities === nothing || length(qualities) == length(sequence) || throw(ArgumentError("qualities and sequence lengths differ"))
    mask, reverse_mask, reverse_shift = _parameters(k); split_mask = ~(UInt64(3) << (k - 1))
    skipped = falses(length(sequence))
    if qualities !== nothing && !allequal(qualities)
        middle = k ÷ 2
        @inbounds for i in k:length(qualities)
            middle_position = i + middle - k + 1
            if Int(qualities[middle_position]) - 33 < minimum_bq
                skipped[i] = true
            end
        end
    end
    result = UInt64[]; forward = UInt64(0); reverse = UInt64(0)
    @inbounds for i in eachindex(sequence)
        nucleotide = encode_base(sequence[i]); complement = UInt64(3) - nucleotide
        forward = ((forward << 2) | nucleotide) & mask
        reverse = ((reverse >> 2) & reverse_mask) | (complement << reverse_shift)
        if i >= k
            split_forward = forward & split_mask; split_reverse = reverse & split_mask
            split_forward == split_reverse && continue
            skipped[i] && continue
            canonical_forward = split_forward < split_reverse
            canonical = canonical_forward ? forward : reverse
            push!(result, canonical | (UInt64(canonical_forward) << 63))
        end
    end
    return result
end

function blockmer_kmers(sequence::AbstractVector{UInt8}, qualities::Union{Nothing,AbstractVector{UInt8}}, k::Integer, l::Integer, minimum_bq::Integer)
    full_length = k + l; length(sequence) < full_length && return Blockmer[]
    k <= 31 && full_length <= 31 || throw(ArgumentError("k and k+l must be <= 31"))
    qualities === nothing || length(qualities) == length(sequence) || throw(ArgumentError("qualities and sequence lengths differ"))
    mask, _, reverse_shift = _parameters(k)
    result = Blockmer[]; forward = UInt64(0); reverse = UInt64(0)
    @inbounds for i in eachindex(sequence)
        nucleotide = encode_base(sequence[i]); complement = UInt64(3) - nucleotide
        forward = ((forward << 2) | nucleotide) & mask
        reverse = (reverse >> 2) | (complement << reverse_shift)
        i < k && continue
        forward == reverse && continue
        if forward < reverse
            i + l <= length(sequence) || continue
            suffix = UInt64(0); low_quality = false
            for position in (i + 1):(i + l)
                suffix = (suffix << 2) | encode_base(sequence[position])
                if qualities !== nothing && Int(qualities[position]) - 33 < minimum_bq
                    low_quality = true; break
                end
            end
            low_quality || push!(result, Blockmer((forward << (2l)) | suffix, true))
        else
            start = i - k + 1; start > l || continue
            suffix = UInt64(0); low_quality = false
            for position in (start - 1):-1:(start - l)
                suffix = (suffix << 2) | (UInt64(3) - encode_base(sequence[position]))
                if qualities !== nothing && Int(qualities[position]) - 33 < minimum_bq
                    low_quality = true; break
                end
            end
            low_quality || push!(result, Blockmer((reverse << (2l)) | suffix, false))
        end
    end
    return result
end
