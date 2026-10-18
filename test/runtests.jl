using Savont
using Test
using CodecZlib

include("rust_oracle.jl")

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
    @test homopolymer_decompress(UInt8['A'], UInt8[1, 2]) == UInt8['A']
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
        @test records[2].sequence == bytes("CGTT")
        @test records[2].qualities == bytes("DCBA")
        @test records[2].reverse_complemented
        @test records[2].sample_index == 2
        compressed = only(read_records(fastq_gz))
        @test compressed.sequence == bytes("ACGT")
        @test compressed.qualities == bytes("IIII")
    end
end
