using Savont
using Test
using CodecZlib

bytes(value) = Vector{UInt8}(codeunits(value))

@testset "configuration" begin
    config = validate(SavontConfig())
    @test config.min_read_length == 1100
    @test config.kmer_size == 17
    @test apply_preset(config, :hifi).min_cluster_size == 4
    @test apply_preset(config, :rrna_operon).max_read_length == 5000
    @test_throws ArgumentError validate(SavontConfig(kmer_size=16))
end

@testset "sequence utilities" begin
    @test reverse_complement("AaCGTn?") == "NNACGTT"
    compressed, lengths = homopolymer_compress(vcat(fill(UInt8('A'), 256), UInt8['C']))
    @test compressed == UInt8['A', 'A', 'C']
    @test lengths == UInt8[255, 1, 1]
    @test homopolymer_decompress(compressed, lengths) == vcat(fill(UInt8('A'), 256), UInt8['C'])
    @test (@test_logs (:warn, r"mismatch") homopolymer_decompress(UInt8['A'], UInt8[1, 2])) == UInt8['A']
    seq, qual, lens = homopolymer_compress_with_quality(bytes("AAAC"), UInt8[70, 60, 65, 50])
    @test (seq, qual, lens) == (bytes("AC"), UInt8[60, 50], UInt8[3, 1])
    @test quality_pool(UInt8[9, 8, 7, 6, 5, 4, 3]) == UInt8[9, 8, 7, 5, 4, 4, 3]
    @test estimate_sequence_identity(UInt8[]) === nothing
    @test estimate_sequence_identity(UInt8[73]) ≈ 99.99
end

@testset "hash and k-mer primitives" begin
    values = UInt64[0, 1, 3, 0x123456789abcdef0, typemax(UInt64)]
    @test all(rev_hash64(mm_hash64(value)) == value for value in values)
    @test decode_kmer(UInt64(0x1b), 4) == "ACGT"
    @test encode_base(UInt8('N')) == 0
    @test fmh_seeds(bytes("ACGTACGTACGT"), 1, 5) |> length == 8
    kmers, positions = fmh_seeds_positions(bytes("ACGTACGTACGT"), 1, 5)
    @test positions == UInt64.(0:7)
    @test all(mm_hash64.(kmers) .== fmh_seeds(bytes("ACGTACGTACGT"), 1, 5))
    @test split_kmer(UInt64(0x1b), 3) == (retrieve_masked_kmer(UInt64(0x1b), 3), UInt64(2))
    constant_qualities = fill(UInt8(33), 12)
    @test split_kmer_mid(bytes("ACGTACGTACGT"), constant_qualities, 5, 25) == split_kmer_mid(bytes("ACGTACGTACGT"), nothing, 5, 25)
    @test !isempty(blockmer_kmers(bytes("AACGTACGTT"), nothing, 5, 2, 25))
    mini, pos = minimizer_seeds_positions(bytes("ACGTACGTACGT"), 3, 5)
    @test length(mini) == length(pos)
    @test !isempty(mini)
end

@testset "FASTA/FASTQ I/O and cutadapt rc" begin
    mktempdir() do directory
        fasta = joinpath(directory, "reads.fa")
        fastq = joinpath(directory, "reads.fastq")
        fastq_gz = joinpath(directory, "compressed.fastq.gz")
        write(fasta, ">read1 ordinary\nACGTN\n")
        write(fastq, "@read2 cutadapt rc\nAACG\n+\nABCD\n")
        open(fastq_gz, "w") do raw
            stream = GzipCompressorStream(raw)
            write(stream, "@read3\nACGT\n+\nIIII\n")
            close(stream)
        end
        records = read_records([fasta, fastq])
        @test length(records) == 2
        @test records[1].sequence == bytes("ACGTN")
        @test records[1].sample_index == 1
        @test (records[1].identifier, records[1].description) == ("read1", "ordinary")
        @test (records[2].identifier, records[2].description) == ("read2", "cutadapt rc")
        @test records[2].sequence == bytes("CGTT")
        @test records[2].qualities == bytes("DCBA")
        @test records[2].reverse_complemented
        @test records[2].sample_index == 2
        compressed = only(read_records(fastq_gz))
        @test compressed.sequence == bytes("ACGT")
        @test compressed.qualities == bytes("IIII")
    end
end

@testset "statistics" begin
    # Values from R: 1 - pbinom(k, n, p) and fisher.test(matrix(...))$p.value.
    @test binomial_test(10, 2, 0.025) ≈ 0.00164317002288 rtol = 1e-10
    @test binomial_test(5, 5, 0.05) == 0.0
    @test fishers_exact((1, 9, 11, 3)) ≈ 0.00275945618522 rtol = 1e-10
    @test fishers_exact((3, 1, 1, 3)) ≈ 0.4857142857 rtol = 1e-9
    @test fishers_exact((0, 0, 4, 5)) == 1.0
    @test_throws ArgumentError fishers_exact((1, 2, 3))
end

@testset "marker inference" begin
    config = SavontConfig(kmer_size=5)
    masked = retrieve_masked_kmer(UInt64(0b0001101100), 5)
    allele(base) = masked | (UInt64(base) << 4)
    counts = KmerCount[(allele(0), (UInt32(40), UInt32(42))), (allele(2), (UInt32(30), UInt32(28))),
                       (allele(3), (UInt32(1), UInt32(1)))]
    info = find_snpmers(counts, config)
    @test info.high_freq_threshold == 100.0
    @test isempty(info.high_freq_kmers)
    snpmer = only(info.snpmer_info)
    @test snpmer.mid_bases == (0x00, 0x02)
    @test snpmer.counts == (UInt32(82), UInt32(58))
    @test snpmer_kmers(snpmer) == (allele(0), allele(2))
    # A minor allele compatible with sequencing error is rejected by the binomial test.
    @test isempty(find_snpmers(KmerCount[(allele(0), (UInt32(500), UInt32(500))), (allele(1), (UInt32(5), UInt32(5)))], config).snpmer_info)
    # Strand bias: the minor allele is seen on one strand only.
    @test isempty(find_snpmers(KmerCount[(allele(0), (UInt32(50), UInt32(50))), (allele(1), (UInt32(60), UInt32(0)))], config).snpmer_info)
    @test isempty(find_snpmers(counts, SavontConfig(kmer_size=5, no_snpmers=true)).snpmer_info)
    @test_throws ArgumentError count_split_kmers(String[], SavontConfig(bloom_filter_size=1.0))
end

# Independent check of the corrected minimizers: every value is the hash of the
# canonical k-mer at its position, and every window minimum is reported inside
# that window.
function check_minimizers(sequence, w, k)
    n_kmers = length(sequence) - k + 1
    values, positions = minimizer_seeds_positions(sequence, w, k)
    n_kmers < w && return isempty(values)
    canonical(position) = begin
        forward = UInt64(0); reverse = UInt64(0)
        for base in @view sequence[position + 1:position + k]
            forward = (forward << 2) | encode_base(base)
        end
        for base in Iterators.reverse(@view sequence[position + 1:position + k])
            reverse = (reverse << 2) | (UInt64(3) - encode_base(base))
        end
        min(forward, reverse)
    end
    hashes = [mm_hash64(canonical(position)) for position in 0:n_kmers - 1]
    all(values .== hashes[Int.(positions) .+ 1]) || return false
    for start in 0:n_kmers - w
        window_minimum = minimum(@view hashes[start + 1:start + w])
        any(i -> values[i] == window_minimum && start <= positions[i] < start + w, eachindex(values)) || return false
    end
    return true
end

@testset "corrected minimizers" begin
    sequences = [record.sequence for record in read_records(joinpath(@__DIR__, "data", "primitives.fasta"))]
    amplicons = read_records(joinpath(@__DIR__, "data", "amplicons_sample1.fastq.gz"))
    append!(sequences, [record.sequence for record in amplicons[1:10]])
    for (w, k) in ((1, 5), (3, 5), (10, 17), (28, 31), (12, 21))
        @test all(sequence -> check_minimizers(sequence, w, k), sequences)
    end
    # The upstream-compatible mode differs only in the first window.
    sequence = amplicons[1].sequence
    fixed_values, fixed_positions = minimizer_seeds_positions(sequence, 28, 31)
    quirk_values, quirk_positions = minimizer_seeds_positions(sequence, 28, 31; upstream_quirks=true)
    @test quirk_values[1] != fixed_values[1]
    @test quirk_values[end-5:end] == fixed_values[end-5:end]
    @test quirk_positions[end-5:end] == fixed_positions[end-5:end]
end

include("rust_oracle.jl")
