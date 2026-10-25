# Differential comparison with the upstream savont 0.8.0 code. When
# SAVONT_RUST_ORACLE points to the compiled oracle (oracle/rust_upstream), its
# output is used directly; otherwise the stored outputs in test/reference are.

using Savont
using Test
using CodecZlib
using Logging

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

optional_list(T, text) = text == "-" ? nothing : parse_list(T, text)

function compare_clusters(text, inputs, config)
    rows = oracle_rows(text)
    paths = joinpath.(DATA_DIRECTORY, inputs)
    info, blockmers = discover_markers(paths, config)
    reads, low_polymorphism = build_twin_reads(paths, info, blockmers, config)
    expected_reads = [fields[2:end] for fields in rows if fields[1] == "read"]
    @test length(reads) == length(expected_reads)
    @test low_polymorphism == (only(only(f[2:end] for f in rows if f[1] == "low_polymorphism")) == "1")
    for (read, f) in zip(reads, expected_reads)
        @test (read.id, read.base_id, read.file_index - 1, read.base_length) ==
              (f[2], f[3], parse(Int, f[4]), parse(Int, f[5]))
        if f[6] == "-"
            @test read.estimated_identity === nothing
        else
            @test read.estimated_identity ≈ parse(Float64, f[6]) rtol = 1e-13
        end
        @test read.sequence == codeunits(f[7])
        @test read.binned_qualities == optional_list(UInt8, f[8])
        @test read.minimizer_positions == parse_list(UInt32, f[9])
        @test read.minimizer_kmers == parse_list(UInt64, f[10])
        @test read.snpmer_positions == parse_list(UInt32, f[11])
        @test read.snpmer_kmers == parse_list(UInt64, f[12])
        @test read.blockmer_positions == parse_list(UInt32, f[13])
        @test read.blockmer_canonical == (parse_list(UInt8, f[14]) .== 1)
        @test read.lsh_signatures == [x == "-" ? nothing : parse(UInt64, x) for x in split(f[15], ',')]
    end
    expected_clusters = [parse_list(Int, fields[2]) .+ 1 for fields in rows if fields[1] == "cluster"]
    clusters = cluster_reads_by_kmers(reads, config)
    @test clusters == expected_clusters
    @test !isempty(expected_clusters)
    low = config.low_polymorphism || low_polymorphism
    if !low
        expected_greedy = [Vector{Int}[] for _ in clusters]
        for fields in rows
            fields[1] == "greedy" && push!(expected_greedy[parse(Int, fields[2]) + 1], parse_list(Int, fields[3]) .+ 1)
        end
        @test greedy_snpmer_clusters(reads, clusters, config) == expected_greedy
    end
    expected_final = [parse_list(Int, fields[2]) .+ 1 for fields in rows if fields[1] == "final"]
    @test cluster_reads_by_snpmers(reads, clusters, config; low_polymorphism=low) == expected_final
    @test !isempty(expected_final)
end

# Per-position pileup summary printed by the oracle `polish` command.
function pileup_summary(pileup)
    join([begin
        matches = count(e -> e.kind == PILEUP_BASE && e.base == p.ref_base, p.entries)
        deletions = count(e -> e.kind == PILEUP_DELETION, p.entries)
        insertions = count(e -> e.kind == PILEUP_INSERTION, p.entries)
        quality = sum((Int(e.quality) for e in p.entries if e.kind == PILEUP_BASE); init=0)
        inserted = sum((Int(e.quality) + (e.inserted == 2 ? Int(e.quality2) : 0) for e in p.entries if e.kind == PILEUP_INSERTION); init=0)
        "$(length(p.entries)):$matches:$deletions:$insertions:$quality:$inserted:$(p.ref_hp_length)"
    end for p in pileup], ',')
end

# Initial consensus (stage 4a), pileups and polishing (stage 4b), merging (stage 5)
# and chimera detection (stage 6).
function compare_consensus(consensus_text, polish_text, merge_text, chimera_text, inputs, config)
    paths = joinpath.(DATA_DIRECTORY, inputs)
    info, blockmers = discover_markers(paths, config)
    reads, low_polymorphism = build_twin_reads(paths, info, blockmers, config)
    clusters = cluster_reads_by_snpmers(reads, cluster_reads_by_kmers(reads, config), config;
                                        low_polymorphism=config.low_polymorphism || low_polymorphism)
    expected = [(parse(Int, f[2]) + 1, parse(Int, f[3]), parse_list(Int, f[4]) .+ 1, f[5]) for f in oracle_rows(consensus_text)]
    consensuses = initial_consensus(reads, clusters, config)
    @test [(c.id, c.depth, c.cluster, String(copy(c.sequence))) for c in consensuses] == expected
    @test !isempty(expected)

    rows = oracle_rows(polish_text)
    by_id = Dict(c.id => c for c in consensuses)
    for f in rows
        f[1] == "pileup" && @test pileup_summary(consensus_pileup(reads, by_id[parse(Int, f[2]) + 1], config)) == f[3]
    end
    kept, low_quality, rates = polish_consensuses!(reads, consensuses, config)
    @test rates == Dict(parse(UInt8, f[2]) => parse(Float64, f[3]) for f in rows if f[1] == "rate")
    low_ids = Set(c.id for c in low_quality)
    polished = sort(vcat(kept, low_quality); by=c -> (-c.depth, c.id))
    @test [(c.id, c.depth, c.id in low_ids, Int.(c.hp_lengths), String(copy(c.sequence)), c.low_quality_positions,
            String(copy(c.decompressed))) for c in polished] ==
          [(parse(Int, f[2]) + 1, parse(Int, f[3]), f[4] == "1", parse_list(Int, f[5]), f[6], parse_list(Int, f[7]) .+ 1, f[8])
           for f in rows if f[1] == "polished"]

    merged = merge_similar_consensuses(kept, config)
    @test [(c.id, c.depth, c.cluster, String(copy(c.decompressed))) for c in merged] ==
          [(parse(Int, f[2]) + 1, parse(Int, f[3]), parse_list(Int, f[4]) .+ 1, f[5]) for f in oracle_rows(merge_text)]

    rows = oracle_rows(chimera_text)
    chimeras = detect_chimeras(merged, config)
    @test [(c.id, string(c.chimera_score)) for c in merged] == [(parse(Int, f[2]) + 1, f[3]) for f in rows if f[1] == "score"]
    @test sort([(merged[c.query].id, merged[c.left_parent].id, merged[c.right_parent].id, c.left_match_length,
                 c.right_match_length, c.query_length, c.coverage_fraction) for c in chimeras]) ==
          sort([(parse(Int, f[2]) + 1, parse(Int, f[3]) + 1, parse(Int, f[4]) + 1, parse(Int, f[5]), parse(Int, f[6]),
                 parse(Int, f[7]), parse(Float64, f[8])) for f in rows if f[1] == "chimera"])
    @test [c.id for c in filter_chimeras(merged, chimeras)] == [parse(Int, f[2]) + 1 for f in rows if f[1] == "kept"]
end

# The whole pipeline: output files, or the early stop with --skip-chimera-detection.
function compare_asv(text, inputs, config)
    rows = split(text, '
'; keepempty=false)
    directory = mktempdir()
    result = savont_asv(joinpath.(DATA_DIRECTORY, inputs), config; output_dir=directory)
    if rows == ["stopped"]
        @test result.asvs === nothing
        @test !isfile(joinpath(directory, "final_asvs.fasta"))
        return
    end
    section(prefix) = [r[length(prefix) + 2:end] for r in rows if startswith(r, prefix * "	")]
    @test readlines(joinpath(directory, "final_asvs.fasta")) == section("fasta")
    @test readlines(joinpath(directory, "feature-table.tsv")) == section("table")
    @test readlines(joinpath(directory, "final_assignments.tsv")) == section("assignment")
    @test !isempty(section("fasta"))
end

function compare_poa(text, input, mode, isa)
    rows = oracle_rows(text)
    records = read_records(joinpath(DATA_DIRECTORY, input))
    sequences = [record.sequence for record in records]
    # As upstream generate_consensus_poa: band base = maximum deviation from the mean length.
    mean_length = sum(length, sequences) ÷ length(sequences)
    deviation = maximum(abs(mean_length - length(sequence)) for sequence in sequences)
    band = mode == "noband" ? nothing : POABand(deviation, 0.1f0)
    graph = POAGraph(); scoring = POAScoring(3, -8, -6)
    for (index, record) in enumerate(records)
        alignment, score = poa_align(graph, record.sequence, scoring; band, isa)
        expected = rows[index]
        @test (expected[1], parse(Int, expected[2])) == ("alignment", index)
        @test score == parse(Int, expected[3])
        @test join(["$(node - 1):$(position - 1)" for (node, position) in alignment], ",") == expected[4]
        poa_add!(graph, alignment, record.sequence, record.qualities)
    end
    @test String(poa_consensus(graph)) == only(fields[2] for fields in rows if fields[1] == "consensus")
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

# Output files of a taxonomy subcommand as `file => lines`, from the oracle text
# or from the files `run(directory)` writes. Rows of equal abundance in the
# single-sample abundance tables are sorted, because upstream aggregates them in
# a hash map and writes ties in arbitrary order.
function taxonomy_files(text::AbstractString)
    files = Dict{String,Vector{String}}()
    for line in split(text, '
'; keepempty=false)
        name, row = split(line, '	'; limit=2)
        push!(get!(Vector{String}, files, name), row)
    end
    for (name, rows) in files
        startswith(rows[1], "abundance	") || continue
        i = 2
        while i <= length(rows)
            j = i
            while j < length(rows) && first(split(rows[j + 1], '	')) == first(split(rows[i], '	'))
                j += 1
            end
            sort!(view(rows, i:j)); i = j + 1
        end
    end
    return files
end

function compare_taxonomy(run, text)
    expected = taxonomy_files(text)
    actual = mktempdir() do directory
        output = joinpath(directory, "out")
        Logging.with_logger(Logging.NullLogger()) do
            run(output)
        end
        taxonomy_files(join(("$name	$line
" for name in sort(readdir(output)) for line in eachline(joinpath(output, name))), ""))
    end
    @test sort(collect(keys(actual))) == sort(collect(keys(expected)))
    for name in keys(expected)
        @test get(actual, name, String[]) == expected[name]
    end
end

source = haskey(ENV, "SAVONT_RUST_ORACLE") ? "live oracle" : "stored oracle output"
@testset "upstream oracle ($source)" begin
    @testset "primitives $input $(p)" for input in PRIMITIVE_INPUTS, p in PRIMITIVE_PARAMETERS
        records = read_records(joinpath(DATA_DIRECTORY, input))
        compare_primitives(oracle_text(o -> primitive_command(o, input, p), primitive_reference_name(input, p)), records, p)
    end
    @testset "partial order alignment $input $mode $isa" for (input, mode, isa) in POA_RUNS
        compare_poa(oracle_text(o -> poa_command(o, input, mode, isa), poa_reference_name(input, mode, isa)), input, mode, isa)
    end
    @testset "statistics" compare_stats(oracle_text(stats_command, STATS_REFERENCE))
    @testset "twin reads and clusters $name" for (name, inputs, arguments, keywords) in CLUSTER_RUNS
        config = SavontConfig(; upstream_quirks=true, keywords...)
        compare_clusters(oracle_text(o -> cluster_command(o, inputs, arguments), cluster_reference_name(name)), inputs, config)
    end
    @testset "consensus, polishing, merging and chimeras $name" for (name, inputs, arguments, keywords) in CLUSTER_RUNS
        config = SavontConfig(; upstream_quirks=true, keywords...)
        compare_consensus(oracle_text(o -> consensus_command(o, inputs, arguments), consensus_reference_name(name)),
                          oracle_text(o -> polish_command(o, inputs, arguments), polish_reference_name(name)),
                          oracle_text(o -> merge_command(o, inputs, arguments), merge_reference_name(name)),
                          oracle_text(o -> chimera_command(o, inputs, arguments), chimera_reference_name(name)), inputs, config)
    end
    @testset "pipeline and output files $name" for (name, inputs, arguments, keywords) in ASV_RUNS
        config = SavontConfig(; upstream_quirks=true, keywords...)
        compare_asv(oracle_text(o -> asv_command(o, inputs, arguments), asv_reference_name(name)), inputs, config)
    end
    @testset "markers $name" for (name, arguments, keywords) in MARKER_RUNS
        config = SavontConfig(; keywords...)
        compare_markers(oracle_text(o -> marker_command(o, arguments), marker_reference_name(name)), config)
    end
    @testset "classify $name" for (name, input, db, arguments, keywords) in CLASSIFY_RUNS
        compare_taxonomy(oracle_text(o -> classify_command(o, input, db, arguments), classify_reference_name(name))) do output
            savont_classify(joinpath(TAXONOMY_DIRECTORY, input), joinpath(TAXONOMY_DIRECTORY, db); output_dir=output, keywords...)
        end
    end
    @testset "sintax $name" for (name, input, db, arguments, keywords) in SINTAX_RUNS
        compare_taxonomy(oracle_text(o -> sintax_command(o, input, db, arguments), sintax_reference_name(name))) do output
            savont_sintax(joinpath(TAXONOMY_DIRECTORY, input), joinpath(TAXONOMY_DIRECTORY, db); output_dir=output, keywords...)
        end
    end
    @testset "export $name" for (name, inputs, arguments, keywords) in EXPORT_RUNS
        compare_taxonomy(oracle_text(o -> export_command(o, inputs, arguments), export_reference_name(name))) do output
            savont_export(joinpath.(TAXONOMY_DIRECTORY, inputs), output; upstream_quirks=true, keywords...)
        end
    end
end
