using Savont
using Test

function parse_csv(text)
    isempty(text) ? UInt64[] : parse.(UInt64, split(text, ','))
end

function rust_oracle_test(executable)
    rows = Dict{String,Vector{Vector{SubString{String}}}}()
    for line in split(read(`$executable`, String), '\n'; keepempty=false)
        fields = split(line, '\t')
        push!(get!(rows, String(fields[1]), Vector{SubString{String}}[]), fields[2:end])
    end
    for row in rows["hash"]
        original, hashed, reversed = parse.(UInt64, row)
        @test mm_hash64(original) == hashed
        @test rev_hash64(hashed) == reversed
    end
    sequence = Vector{UInt8}(codeunits("ACGTACGTACGT"))
    @test fmh_seeds(sequence, 1, 5) == parse_csv(only(rows["fmh"])[1])
    kmers, positions = fmh_seeds_positions(sequence, 1, 5)
    @test kmers == parse_csv(only(rows["fmhp"])[1])
    @test positions == parse_csv(only(rows["fmhp"])[2])
    @test split_kmer_mid(sequence, fill(UInt8(33), length(sequence)), 5, 25) == parse_csv(only(rows["split"])[1])
    expected_blocks = [(parse(UInt64, first(split(item, ':'))), last(split(item, ':')) == "1") for item in split(only(rows["blocks"])[1], ',')]
    @test [(x.kmer, x.is_forward) for x in blockmer_kmers(Vector{UInt8}(codeunits("AACGTACGTT")), nothing, 5, 2, 25)] == expected_blocks
    minimizers, minimizer_positions = minimizer_seeds_positions(sequence, 3, 5)
    @test minimizers == parse_csv(only(rows["minimizer"])[1])
    @test minimizer_positions == parse_csv(only(rows["minimizer"])[2])
    @test reverse_complement("AaCGTn?") == only(rows["revcomp"])[1]
end

if haskey(ENV, "SAVONT_RUST_ORACLE")
    @testset "Rust differential oracle" rust_oracle_test(ENV["SAVONT_RUST_ORACLE"])
end
