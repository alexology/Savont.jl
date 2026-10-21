# Differential comparison with the upstream savont 0.8.0 code. When
# SAVONT_RUST_ORACLE points to the compiled oracle (oracle/rust_upstream), its
# output is used directly; otherwise the stored outputs in test/reference are.

using Savont
using Test
using CodecZlib

include("oracle_cases.jl")

parse_list(T, text) = isempty(text) ? T[] : parse.(T, split(text, ','))

function oracle_text(command, reference)
    oracle = get(ENV, "SAVONT_RUST_ORACLE", "")
    isempty(oracle) || return read(command(oracle), String)
    return String(transcode(GzipDecompressor, read(joinpath(REFERENCE_DIRECTORY, reference))))
end

oracle_rows(text) = [split(line, '\t') for line in split(text, '\n'; keepempty=false)]

function compare_primitives(text, records, p)
    rows = Dict{Tuple{String,Int},Vector{SubString{String}}}()
    for fields in oracle_rows(text)
        rows[(String(fields[1]), parse(Int, fields[2]))] = fields[3:end]
    end
    for (index, record) in enumerate(records)
        sequence = record.sequence; qualities = record.qualities
        @test fmh_seeds(sequence, p.c, p.k) == parse_list(UInt64, rows[("fmh", index)][1])
        kmers, positions = fmh_seeds_positions(sequence, p.c, p.k)
        @test kmers == parse_list(UInt64, rows[("fmhp", index)][1])
        @test positions == parse_list(UInt64, rows[("fmhp", index)][2])
        minimizers, minimizer_positions = minimizer_seeds_positions(sequence, p.w, p.k; upstream_quirks=true)
        @test minimizers == parse_list(UInt64, rows[("minimizer", index)][1])
        @test minimizer_positions == parse_list(UInt64, rows[("minimizer", index)][2])
        @test split_kmer_mid(sequence, qualities, p.k, p.minimum_bq) == parse_list(UInt64, rows[("split", index)][1])
        expected_blocks = [(parse(UInt64, first(split(item, ':'))), last(item) == '1')
                           for item in split(rows[("blocks", index)][1], ','; keepempty=false)]
        @test [(x.kmer, x.is_forward) for x in blockmer_kmers(sequence, qualities, p.k, p.l, p.minimum_bq)] == expected_blocks
        @test reverse_complement(sequence) == codeunits(rows[("revcomp", index)][1])
        compressed, lengths = homopolymer_compress(sequence)
        @test compressed == codeunits(rows[("hpc", index)][1])
        @test lengths == parse_list(UInt8, rows[("hpc", index)][2])
        if qualities !== nothing
            compressed, compressed_qualities, lengths = homopolymer_compress_with_quality(sequence, qualities)
            expected = rows[("hpcq", index)]
            @test (compressed, compressed_qualities, lengths) ==
                  (codeunits(expected[1]), parse_list(UInt8, expected[2]), parse_list(UInt8, expected[3]))
            @test quality_pool(qualities) == parse_list(UInt8, rows[("pool", index)][1])
            @test estimate_sequence_identity(qualities) ≈ parse(Float64, rows[("identity", index)][1]) rtol = 1e-13
        end
    end
end

function compare_markers(text, config)
    rows = oracle_rows(text)
    select(kind) = [fields[2:end] for fields in rows if fields[1] == kind]
    as_counts(selected) = [(parse(UInt64, f[1]), (parse(UInt32, f[2]), parse(UInt32, f[3]))) for f in selected]
    paths = joinpath.(DATA_DIRECTORY, AMPLICON_INPUTS)
    snpmer_counts, blockmer_counts = count_split_kmers(paths, config)
    @test snpmer_counts == as_counts(select("count"))
    @test blockmer_counts == as_counts(select("blockcount"))
    info, blockmers = discover_markers(paths, config)
    @test info.high_freq_threshold == parse(Float64, only(only(select("high_freq_threshold"))))
    @test sort!(collect(info.high_freq_kmers)) == [parse(UInt64, only(f)) for f in select("high_freq")]
    expected_snpmers = [SnpmerInfo(parse(UInt64, f[1]), (parse(UInt8, f[2]), parse(UInt8, f[3])),
                                   (parse(UInt32, f[4]), parse(UInt32, f[5])), UInt8(config.kmer_size))
                        for f in select("snpmer")]
    @test info.snpmer_info == expected_snpmers
    @test !isempty(expected_snpmers)
    expected_blockmers = [(parse(UInt64, f[1]), parse(UInt64, f[2]), parse(UInt64, f[3]),
                           parse(UInt32, f[4]), parse(UInt32, f[5])) for f in select("blockmer")]
    @test [(b.anchor_kmer, b.blockmers[1].kmer, b.blockmers[2].kmer, b.counts...) for b in blockmers] == expected_blockmers
    config.use_blockmers && @test !isempty(expected_blockmers)
end

function compare_stats(text)
    binomial_cases = 0; fisher_cases = 0
    for fields in oracle_rows(text)
        if fields[1] == "binomial"
            n, k = parse.(Int, fields[2:3])
            @test binomial_test(n, k, parse(Float64, fields[4])) ≈ parse(Float64, fields[5]) rtol = 1e-12 atol = 1e-300
            binomial_cases += 1
        else
            table = parse.(Int, fields[2:5])
            @test fishers_exact(table) ≈ parse(Float64, fields[6]) rtol = 1e-12 atol = 1e-300
            fisher_cases += 1
        end
    end
    @test binomial_cases > 1000 && fisher_cases > 1000
end

source = haskey(ENV, "SAVONT_RUST_ORACLE") ? "live oracle" : "stored oracle output"
@testset "upstream oracle ($source)" begin
    @testset "primitives $input $(p)" for input in PRIMITIVE_INPUTS, p in PRIMITIVE_PARAMETERS
        records = read_records(joinpath(DATA_DIRECTORY, input))
        compare_primitives(oracle_text(o -> primitive_command(o, input, p), primitive_reference_name(input, p)), records, p)
    end
    @testset "statistics" compare_stats(oracle_text(stats_command, STATS_REFERENCE))
    @testset "markers $name" for (name, arguments, keywords) in MARKER_RUNS
        config = SavontConfig(; keywords...)
        compare_markers(oracle_text(o -> marker_command(o, arguments), marker_reference_name(name)), config)
    end
end
