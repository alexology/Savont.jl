# Phase 3C, first part: retained reads (upstream `seeding::get_twin_read_syncmer`,
# `kmer_comp::twin_reads_from_snpmers` and `main.rs::get_twin_reads_from_kmer_info`).

const MAX_KMER_COUNT_IN_READ = 500
const QUALITY_SEQ_BIN = 4
const LSH_NUM_TABLES = 20
const LSH_BUCKET_SIZE = 3

"""
    TwinRead

A retained read with its markers (upstream `TwinRead`). Positions are 0-based
k-mer start positions, as upstream.

- `id`: full header line; `base_id`: its first word; `file_index`: 1-based input
  file; `base_length`: read length.
- `sequence`: upper-case bases, with every non-ACGT byte (including `N` and `U`)
  replaced by `A`, as upstream's 2-bit storage.
- `binned_qualities`: for each block of four bases, the minimum Phred+33 byte
  rounded to upstream's 4-bit `QualCompact3` code (`(byte - 32) ÷ 3`, clamped to
  0–15; the represented byte is `3 * code + 33`); `nothing` for FASTA input.
- `estimated_identity`: mean per-base accuracy in percent from the qualities;
  `nothing` without qualities or when all qualities are equal.
- `minimizer_kmers`/`minimizer_positions`: open syncmers (see [`twin_read`](@ref)).
  As upstream, filtering of high-frequency k-mers removes entries from the
  positions only: the k-mer vectors keep every syncmer and SNPmer and are what
  primary clustering uses.
- `snpmer_kmers`/`snpmer_positions`: occurrences of the inferred SNPmers, with the
  same position-only filtering.
- `blockmer_positions`/`blockmer_canonical`: occurrences of inferred blockmers and
  whether the forward (`true`) or reverse orientation matched.
- `lsh_signatures`: one bucket signature per LSH table, computed from
  `minimizer_kmers`; `nothing` with fewer than three minimizers.
"""
mutable struct TwinRead
    id::String
    base_id::String
    file_index::Int
    base_length::Int
    sequence::Vector{UInt8}
    binned_qualities::Union{Nothing,Vector{UInt8}}
    estimated_identity::Union{Nothing,Float64}
    k::UInt8
    l::UInt8
    minimizer_positions::Vector{UInt32}
    minimizer_kmers::Vector{UInt64}
    snpmer_positions::Vector{UInt32}
    snpmer_kmers::Vector{UInt64}
    blockmer_positions::Vector{UInt32}
    blockmer_canonical::Vector{Bool}
    lsh_signatures::Vector{Union{Nothing,UInt64}}
end

_first_word(text::AbstractString) = (words = split(text); isempty(words) ? String(text) : String(first(words)))

function _dna_sequence(sequence::AbstractVector{UInt8})
    return map(sequence) do base
        upper = UInt8(uppercase(Char(base)))
        upper in (UInt8('A'), UInt8('C'), UInt8('G'), UInt8('T')) ? upper : UInt8('A')
    end
end

function _binned_qualities(qualities::AbstractVector{UInt8})
    binned = UInt8[]
    for start in 1:QUALITY_SEQ_BIN:length(qualities)
        lowest = minimum(@view qualities[start:min(start + QUALITY_SEQ_BIN - 1, end)])
        push!(binned, UInt8(clamp((Int(lowest) - 32) ÷ 3, 0, 15)))
    end
    return binned
end

"""
    twin_read(sequence, qualities, id, config, snpmer_set, blockmer_set; file_index=1) -> Union{TwinRead,Nothing}

Build the markers of one read (upstream `get_twin_read_syncmer`); `nothing` when
the read is shorter than `k`. Minimizers are open syncmers: a canonical k-mer is
selected when, among its `k - s + 1` canonical s-mers (`s = k - sampling_rate + 1`),
the middle one has the strictly smallest `mm_hash64`. SNPmers are canonical k-mers
in `snpmer_set` whose middle base quality exceeds `minimum_base_quality` (or whose
read has constant qualities); a SNPmer context seen more than once in the read is
dropped. Blockmers are `k + l`-mers in `blockmer_set` whose `l` suffix bases all
have quality above `minimum_base_quality`.

With `config.upstream_quirks = false` (default), syncmers are selected along the
whole read. Upstream initializes its rolling s-mer with only the first `s - 1`
bases, so the first `k - s` k-mers of a read are never examined and the next
`s - 1` are judged on hashes of non-contiguous bases (deviation D2 in
`docs/deviations.md`); `upstream_quirks = true` reproduces that.
"""
function twin_read(sequence::AbstractVector{UInt8}, qualities::Union{Nothing,AbstractVector{UInt8}},
                   id::AbstractString, config::SavontConfig,
                   snpmer_set::AbstractSet{UInt64}, blockmer_set::AbstractSet{UInt64}; file_index::Integer=1)
    k = config.kmer_size; c = config.sampling_rate; l = config.blockmer_length
    minimum_bq = Int(config.minimum_base_quality); quirks = config.upstream_quirks
    n = length(sequence)
    n < k && return nothing
    qualities === nothing || length(qualities) == n ||
        throw(ArgumentError("qualities and sequence lengths differ for read $id"))
    s = k - c + 1
    s >= 1 || throw(ArgumentError("sampling_rate must not exceed kmer_size"))
    blockmer_k = k + l
    mask, reverse_mask, reverse_shift = _parameters(k)
    s_mask, s_reverse_mask, s_reverse_shift = _parameters(s)
    blockmer_mask, blockmer_reverse_mask, blockmer_reverse_shift = _parameters(blockmer_k)
    split_mask = ~(UInt64(3) << (k - 1))
    middle_offset = k ÷ 2
    window = k - s + 1
    middle_index = (k - s) ÷ 2 + 1
    equal_qualities = qualities !== nothing && allequal(qualities)

    forward = UInt64(0); reverse = UInt64(0)
    s_forward = UInt64(0); s_reverse = UInt64(0)
    block_forward = UInt64(0); block_reverse = UInt64(0)
    s_hashes = UInt64[]
    minimizer_positions = UInt32[]; minimizer_kmers = UInt64[]
    snpmer_positions = UInt32[]; snpmer_kmers = UInt64[]
    blockmer_positions = UInt32[]; blockmer_canonical = Bool[]
    context_counts = Dict{UInt64,Int}()

    # Bases are 0-based (`i`) in the loops below, as upstream.
    for i in 0:(k - 2)
        nucleotide = encode_base(sequence[i + 1]); complement = UInt64(3) - nucleotide
        forward = (forward << 2) | nucleotide
        reverse = (reverse >> 2) | (complement << reverse_shift)
        if quirks
            if i < s - 1
                s_forward = (s_forward << 2) | nucleotide
                s_reverse = (s_reverse >> 2) | (complement << s_reverse_shift)
            end
        else
            s_forward = ((s_forward << 2) | nucleotide) & s_mask
            s_reverse = ((s_reverse >> 2) & s_reverse_mask) | (complement << s_reverse_shift)
            if i >= s - 1
                push!(s_hashes, mm_hash64(min(s_forward, s_reverse)))
                length(s_hashes) > window && popfirst!(s_hashes)
            end
        end
    end
    for i in 0:min(blockmer_k - 2, n - 1)
        nucleotide = encode_base(sequence[i + 1])
        block_forward = (block_forward << 2) | nucleotide
        block_reverse = (block_reverse >> 2) | ((UInt64(3) - nucleotide) << blockmer_reverse_shift)
    end

    for i in (k - 1):(n - 1)
        nucleotide = encode_base(sequence[i + 1]); complement = UInt64(3) - nucleotide
        forward = ((forward << 2) | nucleotide) & mask
        reverse = ((reverse >> 2) & reverse_mask) | (complement << reverse_shift)
        canonical = (forward & split_mask) < (reverse & split_mask) ? forward : reverse

        s_forward = ((s_forward << 2) | nucleotide) & s_mask
        s_reverse = ((s_reverse >> 2) & s_reverse_mask) | (complement << s_reverse_shift)
        push!(s_hashes, mm_hash64(min(s_forward, s_reverse)))
        length(s_hashes) > window && popfirst!(s_hashes)

        if i >= blockmer_k - 1
            block_forward = ((block_forward << 2) | nucleotide) & blockmer_mask
            block_reverse = ((block_reverse >> 2) & blockmer_reverse_mask) | (complement << blockmer_reverse_shift)
            forward_found = block_forward in blockmer_set
            if forward_found || block_reverse in blockmer_set
                start = i + 1 - blockmer_k
                suffix_good = true
                if qualities !== nothing && !equal_qualities
                    for j in 0:(l - 1)
                        position = start + k + j
                        if position < length(qualities) && Int(qualities[position + 1]) - 33 <= minimum_bq
                            suffix_good = false
                            break
                        end
                    end
                end
                if suffix_good || equal_qualities
                    push!(blockmer_positions, UInt32(start)); push!(blockmer_canonical, forward_found)
                end
            end
        end

        if canonical in snpmer_set
            middle_quality = qualities === nothing ? 60 : Int(qualities[i + 1 + middle_offset - k + 1]) - 33
            if middle_quality > minimum_bq || equal_qualities
                push!(snpmer_positions, UInt32(i + 1 - k)); push!(snpmer_kmers, canonical)
            end
            context = canonical & split_mask
            context_counts[context] = get(context_counts, context, 0) + 1
        end

        if length(s_hashes) == window
            middle_hash = s_hashes[middle_index]
            if all(j -> j == middle_index || s_hashes[j] > middle_hash, 1:window)
                push!(minimizer_positions, UInt32(i + 1 - k)); push!(minimizer_kmers, canonical)
            end
        end
    end

    unique_context = [context_counts[kmer & split_mask] == 1 for kmer in snpmer_kmers]
    identity = qualities === nothing || equal_qualities ? nothing : estimate_sequence_identity(qualities)
    return TwinRead(String(id), _first_word(id), file_index, n, _dna_sequence(sequence),
                    qualities === nothing ? nothing : _binned_qualities(qualities), identity,
                    UInt8(k), UInt8(l), minimizer_positions, minimizer_kmers,
                    snpmer_positions[unique_context], snpmer_kmers[unique_context],
                    blockmer_positions, blockmer_canonical, Union{Nothing,UInt64}[])
end

# rustc's FxHasher64 (fxhash 0.2.1) with one `write_u64` per word.
@inline _fx_hash_word(hash::UInt64, word::UInt64) = (bitrotate(hash, 5) ⊻ word) * 0x517cc1b727220a95

"""
    lsh_signatures(kmers) -> Vector{Union{Nothing,UInt64}}

MinHash bucket signatures (upstream `TwinRead::compute_lsh_signatures`): for each
of 20 tables, the k-mers are ordered by an FxHash seeded with the table index and
the three smallest are combined as `xor(kmer_i * i)`, `i = 1, 2, 3`.
"""
function lsh_signatures(kmers::AbstractVector{UInt64})
    length(kmers) < LSH_BUCKET_SIZE && return Union{Nothing,UInt64}[nothing for _ in 1:LSH_NUM_TABLES]
    signatures = Vector{Union{Nothing,UInt64}}(undef, LSH_NUM_TABLES)
    for table in 0:(LSH_NUM_TABLES - 1)
        seeded = _fx_hash_word(UInt64(0), UInt64(table))
        order = sortperm([_fx_hash_word(seeded, kmer) for kmer in kmers])  # stable, as sort_by_key
        signature = UInt64(0)
        for i in 1:LSH_BUCKET_SIZE
            signature ⊻= kmers[order[i]] * UInt64(i)
        end
        signatures[table + 1] = signature
    end
    return signatures
end

"""
    build_twin_reads(paths, info::KmerGlobalInfo, blockmers, config) -> (reads, low_polymorphism)

Build the retained reads of `paths` (upstream `twin_reads_from_snpmers` followed
by the sorting in `get_twin_reads_from_kmer_info`):

1. reads outside `min_read_length:max_read_length` are skipped (a warning is
   emitted when this removes more than half of a file); unlike k-mer counting,
   Cutadapt `rc` reads are used as written;
2. markers are extracted with [`twin_read`](@ref); syncmers occurring more than
   500 times in the read or in the high-frequency set, and SNPmers in that set,
   are removed from the positions only (see [`TwinRead`](@ref));
3. reads keeping fewer than `length ÷ sampling_rate ÷ 20` syncmer positions are
   dropped as repetitive;
4. reads are sorted by header, those with an estimated identity below
   `quality_value_cutoff` are dropped, and the rest are stably sorted by
   decreasing estimated identity (reads without one count as 100).

`low_polymorphism` is `true` when more than 75% of the reads have no SNPmer, the
condition under which upstream enables its low-polymorphism mode.
"""
function build_twin_reads(paths::AbstractVector{<:AbstractString}, info::KmerGlobalInfo,
                          blockmers::AbstractVector{BlockmerInfo}, config::SavontConfig)
    validate(config)
    snpmer_set = Set{UInt64}()
    for snpmer in info.snpmer_info
        union!(snpmer_set, snpmer_kmers(snpmer))
    end
    blockmer_set = Set{UInt64}(b.kmer for info in blockmers for b in info.blockmers)
    high_freq = info.high_freq_kmers
    reads = TwinRead[]
    seen = zeros(Int, length(paths)); length_removed = zeros(Int, length(paths))
    foreach_record(paths; cutadapt_rc=false) do record
        sample = record.sample_index
        seen[sample] += 1
        n = length(record.sequence)
        if n < config.min_read_length || n > config.max_read_length
            length_removed[sample] += 1
            return
        end
        read = twin_read(record.sequence, record.qualities, record.header, config, snpmer_set, blockmer_set;
                         file_index=sample)
        read === nothing && return
        occurrences = Dict{UInt64,Int}()
        for kmer in read.minimizer_kmers
            occurrences[kmer] = get(occurrences, kmer, 0) + 1
        end
        solid = [occurrences[kmer] <= MAX_KMER_COUNT_IN_READ && !(kmer in high_freq) for kmer in read.minimizer_kmers]
        count(solid) < n ÷ config.sampling_rate ÷ 20 && return
        read.minimizer_positions = read.minimizer_positions[solid]
        read.snpmer_positions = read.snpmer_positions[[!(kmer in high_freq) for kmer in read.snpmer_kmers]]
        read.lsh_signatures = lsh_signatures(read.minimizer_kmers)
        push!(reads, read)
    end
    for (sample, path) in enumerate(paths)
        if length_removed[sample] > seen[sample] ÷ 2
            @warn "More than 50% of reads were removed in $path by length filtering " *
                  "(min: $(config.min_read_length), max: $(config.max_read_length))"
        end
    end
    sort!(reads; by=read -> read.id)
    below = count(read -> read.estimated_identity !== nothing &&
                          read.estimated_identity < config.quality_value_cutoff, reads)
    if below / length(reads) > 0.5
        @warn "More than 50% of reads are below the quality threshold of $(config.quality_value_cutoff)%"
    end
    filter!(read -> read.estimated_identity === nothing || read.estimated_identity >= config.quality_value_cutoff, reads)
    sort!(reads; by=read -> -something(read.estimated_identity, 100.0))
    fraction_without = count(read -> isempty(read.snpmer_positions), reads) / length(reads)
    if fraction_without > 0.10
        @warn "High fraction of reads without SNPmers: $(round(100 * fraction_without; digits=2))%. " *
              "This may indicate low polymorphism in the sample."
    end
    return reads, fraction_without > 0.75
end
