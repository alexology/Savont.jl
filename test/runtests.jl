using Savont
using Test
using CodecZlib
using Logging
using Tar

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
    @test homopolymer_compress(bytes("GAAT"), false) == (bytes("GAAT"), UInt8[1, 1, 1, 1])
    @test homopolymer_compress(bytes("GAAT"), false; upstream_quirks=true) == (bytes("GAATG"), UInt8[1, 1, 1, 1, 1])  # D4
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
    # D21 (SV17): extreme counts, against R `pbinom(k, n, p, lower.tail = FALSE)`.
    @test binomial_test(1_000_000_000, 25_000_000, 0.025) ≈ 0.499946803434366 rtol = 1e-5
    @test binomial_test(100_000_000, 5_000_000, 0.05) ≈ 0.499881 rtol = 1e-5
    @test binomial_test(1_000_000_000, 25_000_000, 0.025; upstream_quirks=true) == 1.1876349336853742   # statrs
    @test binomial_test(100_000, 2_500, 0.025) == binomial_test(100_000, 2_500, 0.025; upstream_quirks=true)
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

@testset "syncmers along the whole read" begin
    no_markers = Set{UInt64}()
    syncmers(sequence; quirks=false) = begin
        read = twin_read(sequence, nothing, "read", SavontConfig(upstream_quirks=quirks), no_markers, no_markers)
        Dict(Int(p) => kmer for (p, kmer) in zip(read.minimizer_positions, read.minimizer_kmers))
    end
    # A syncmer depends only on the k-mer sequence: trimming the read start only
    # shifts the syncmers of the remaining sequence.
    for record in read_records(joinpath(@__DIR__, "data", "amplicons_sample1.fastq.gz"))[1:30], shift in (1, 3, 8)
        full = syncmers(record.sequence); trimmed = syncmers(record.sequence[shift + 1:end])
        @test trimmed == Dict(p - shift => kmer for (p, kmer) in full if p >= shift)
    end
    sequence = read_records(joinpath(@__DIR__, "data", "amplicons_sample1.fastq.gz"))[1].sequence
    @test minimum(keys(syncmers(sequence; quirks=true))) >= 10  # upstream: k - s = 10
    @test syncmers(sequence) != syncmers(sequence; quirks=true)
    @test twin_read(sequence[1:16], nothing, "short", SavontConfig(), no_markers, no_markers) === nothing
end

@testset "low-polymorphism primary clusters (D3)" begin
    paths = joinpath.(@__DIR__, "data", ["amplicons_sample1.fastq.gz", "amplicons_sample2.fastq.gz"])
    final_sizes = map((true, false)) do quirks
        config = SavontConfig(upstream_quirks=quirks)
        info, blockmers = discover_markers(paths, config)
        reads, low = @test_logs (:warn, r"without SNPmers") match_mode=:any build_twin_reads(paths, info, blockmers, config)
        @test !low
        clusters = cluster_reads_by_kmers(reads, config)
        final = cluster_reads_by_snpmers(reads, clusters, config)
        # The unrelated amplicon has no SNPmers and forms the 45-read primary cluster.
        unrelated = only(filter(c -> all(i -> isempty(reads[i].snpmer_positions), c), clusters))
        @test length(unrelated) == 45
        @test (unrelated in final) == !quirks
        length.(final)
    end
    @test final_sizes == ([112, 75, 25], [112, 75, 45, 25])
end

@testset "k-mer deduplication of consensuses (D5)" begin
    sequence = read_records(joinpath(@__DIR__, "data", "noisy_1pct.fastq.gz"))[1].sequence
    consensus(id, depth, s) = ConsensusSequence(id, copy(s), ones(UInt8, length(s)), depth, Int[])
    # The second is the first trimmed at both ends, with less than half its depth.
    consensuses = [consensus(1, 100, sequence), consensus(2, 40, sequence[11:end-10]), consensus(3, 60, sequence[11:end-10]),
                   consensus(4, 10, sequence[1:90])]
    @test [c.id for c in remove_similar_consensuses(consensuses)] == [1, 3]
    # Upstream: kept consensuses are zeroed before the later comparisons, so nothing
    # but the consensus shorter than 100 bases is removed.
    @test [c.id for c in remove_similar_consensuses(consensuses; upstream_quirks=true)] == [1, 2, 3]
end

@testset "stage 5 with only short consensuses (D7)" begin
    short = ConsensusSequence(1, Vector{UInt8}(codeunits(repeat("ACGT", 20))), ones(UInt8, 80), 12, collect(1:12))
    @test isempty(remove_similar_consensuses([short]))
    @test isempty(merge_similar_consensuses([short], SavontConfig()))
    # Upstream panics when building the minimap2 index of no sequence
    @test_throws ErrorException merge_similar_consensuses([short], SavontConfig(upstream_quirks=true))
end

@testset "chimera similarity lookup (D6)" begin
    table = Dict((3, 1) => 0.9)        # keyed (later, earlier) as upstream
    @test Savont._similarity(table, 1, 3, 0.5, false) == 0.9
    @test Savont._similarity(table, 3, 1, 0.5, false) == 0.9
    @test Savont._similarity(table, 1, 3, 0.5, true) == 0.5   # upstream: never found
    @test Savont._as_i64(-Inf) == typemin(Int) && Savont._as_i64(-12.9) == -12 && Savont._as_i64(NaN) == 0
end

@testset "skipping chimera detection completes the run (D8), tagged read headers (D9)" begin
    paths = joinpath.(@__DIR__, "data", ["amplicons_sample1.fastq.gz", "amplicons_sample2.fastq.gz"])
    directory = mktempdir()
    result = savont_asv(paths, SavontConfig(skip_chimera_detection=true); output_dir=directory)
    @test !isempty(result.asvs) && isempty(result.chimeras)
    @test all(isfile(joinpath(directory, f)) for f in ("final_asvs.fasta", "feature-table.tsv", "final_assignments.tsv"))
    @test savont_asv(paths, SavontConfig(skip_chimera_detection=true, upstream_quirks=true)).asvs === nothing

    # D9: headers with tab-separated tags (as written by dorado).
    tagged = deepcopy(result.reads)
    foreach(read -> read.id *= "\tqs:f:20.5\tch:i:7", tagged)
    function assignments(reads, quirks)
        path = joinpath(directory, "assignments.tsv")
        write_final_assignments(path, reads, result.asvs, result.mappings; upstream_quirks=quirks)
        return readlines(path)
    end
    plain = assignments(result.reads, false)
    @test assignments(tagged, false) == plain
    @test count(line -> !occursin("\tunassigned\t", line), plain) > 200
    # Upstream: every tagged read is unassigned, and the header keeps its tags.
    @test all(line -> occursin("\tch:i:7\tunassigned\tNA\t", line), assignments(tagged, true)[2:end])
end

@testset "odd sampling rate (D11), conserving EM depths (D12), duplicate read ids (D13)" begin
    # D11: with an even sampling rate the two strands share no syncmer (D2 corrected,
    # so that read ends do not differ; `twin_read` does not validate the rate).
    rng_sequence = [b"ACGT"[(mm_hash64(UInt64(i)) & 3) + 1] for i in 1:1300]
    syncmers(c, sequence) = Set(twin_read(sequence, nothing, "r", SavontConfig(sampling_rate=c),
                                          Set{UInt64}(), Set{UInt64}()).minimizer_kmers)
    @test syncmers(11, rng_sequence) == syncmers(11, reverse_complement(rng_sequence))
    @test isempty(intersect(syncmers(10, rng_sequence), syncmers(10, reverse_complement(rng_sequence))))
    @test_throws ArgumentError validate(SavontConfig(sampling_rate=10))
    @test validate(SavontConfig(sampling_rate=10, upstream_quirks=true)).sampling_rate == 10

    # D12: depths of ambiguous reads add up to the assigned reads.
    for (classes, n, upstream, corrected) in (([[1, 2] => 1], 2, [1, 1], [1, 0]), ([[1, 2, 3] => 1], 3, [0, 0, 0], [1, 0, 0]),
                                              ([[1, 2] => 31], 2, [16, 16], [16, 15]), ([[1] => 5, [1, 2] => 4, [2] => 2], 2, nothing, nothing))
        assigned = sum(last, classes)
        quirky = Savont._em_depths(classes, n, assigned; upstream_quirks=true)
        fixed = Savont._em_depths(classes, n, assigned)
        upstream === nothing || @test quirky == upstream
        corrected === nothing || @test fixed == corrected
        @test sum(fixed) == assigned
        @test all(abs.(fixed .- quirky) .<= 1)
    end

    # D13: two reads with the same id keep their own assignments.
    paths = joinpath.(@__DIR__, "data", ["amplicons_sample1.fastq.gz", "amplicons_sample2.fastq.gz"])
    result = Logging.with_logger(Logging.NullLogger()) do
        savont_asv(paths, SavontConfig(pooled_samples=true))
    end
    final = Set(c.id for c in result.asvs)
    assigned_to(r) = (i = findfirst(m -> m.asv_id in final, result.mappings[r]); i === nothing ? nothing : result.mappings[r][i].asv_id)
    a = findfirst(r -> assigned_to(r) !== nothing, eachindex(result.reads))
    b = findfirst(r -> assigned_to(r) !== nothing && assigned_to(r) != assigned_to(a) && result.reads[r].file_index == 2, eachindex(result.reads))
    reads = deepcopy(result.reads); reads[b].id = reads[a].id
    rows(quirks) = mktempdir() do directory
        path = joinpath(directory, "assignments.tsv")
        write_final_assignments(path, reads, result.asvs, result.mappings; upstream_quirks=quirks)
        split.(readlines(path)[[a + 1, b + 1]], '\t')
    end
    @test rows(false)[1][2] != rows(false)[2][2]
    @test rows(true)[1][2] == rows(true)[2][2]
end

const TAXONOMY_DATA = joinpath(@__DIR__, "data", "taxonomy")

@testset "taxonomy primitives" begin
    @test Savont._rust_float.([1.0, 0.0, 0.5, 0.42000000000000004, 1.0e-5, 12.5, 1.0e20, 100.0]) ==
          ["1", "0", "0.5", "0.42000000000000004", "0.00001", "12.5", "100000000000000000000", "100"]
    rng = Savont._Xorshift(0)                       # seeded with max(seed, 1)
    @test Savont._next!(rng) == 1082269761
    @test isempty(sintax_kmers(bytes("ACGTACGTACG"))) && length(sintax_kmers(bytes("ACGTACGTACGTA"))) == 2
    @test sintax_kmers(bytes("ACGUACGUACGU")) == sintax_kmers(bytes("acgtacgtacgt"))
    @test sequence_hash(bytes("ACGTTGCAACGG")) == sequence_hash(reverse_complement(bytes("acgttgcaacgg")))
    @test length(sequence_hash(bytes("A"))) == 16
    entry = TaxonomyEntry("7", "S s", "G", "F", "O", "C", "P", "Cl", "K", "sub", "", "")
    levels(identity) = (t = assign_taxonomy(entry, identity, 99.0, 94.5, "asv1", false);
                        [t.species, t.genus, t.family, t.order, t.class, t.phylum, t.superkingdom])
    @test levels(99.0) == ["S s", "G", "F", "O", "C", "P", "K"]
    @test levels(94.5) == ["UNCLASSIFIED", "G", "F", "O", "C", "P", "K"]
    @test levels(86.5)[2:3] == ["UNCLASSIFIED", "F"] && levels(82.0)[3:4] == ["UNCLASSIFIED", "O"]
    @test levels(78.5)[4:5] == ["UNCLASSIFIED", "C"] && levels(75.0)[5:6] == ["UNCLASSIFIED", "P"]
    @test levels(74.9) == fill("UNCLASSIFIED", 7)
    @test assign_taxonomy(entry, 50, 99.0, 94.5, "asv1", true).clade == "UNCLASSIFIED-(asv1)"
    @test assign_taxonomy(entry, 99.5, 99.0, 94.5, "asv1", false).subspecies == "sub"
    @test assign_taxonomy(entry, 95, 99.0, 94.5, "asv1", false).subspecies == ""
end

@testset "reference databases" begin
    emu = Logging.with_logger(Logging.NullLogger()) do
        load_database(joinpath(TAXONOMY_DATA, "emu-1"))
    end
    @test emu.keyword == "emu-1" && length(emu.taxonomy) == 16    # malformed line skipped
    @test database_key(emu, ">1001:emu_db:1 [x]") == "1001"
    silva = Logging.with_logger(() -> load_database(joinpath(TAXONOMY_DATA, "silva-138.2")), Logging.NullLogger())
    @test database_key(silva, "AB000001.1.1050 Bacteria;Bacillota") == "AB000001" && database_key(silva, "  ") === nothing
    @test silva.taxonomy["AB000001"].genus == "Lactobacillus" && silva.taxonomy["AB000001"].species == "Lactobacillus acidophilus"
    short_path = only(e for e in values(silva.taxonomy) if e.order == "Synechococcales")
    @test short_path.family == "" && short_path.genus == "UNKNOWN"
    gg2 = load_database(joinpath(TAXONOMY_DATA, "greengenes2-2024.09"))
    @test any(e.species == "Lactobacillus acidophilus group" for e in values(gg2.taxonomy))
    @test any(e.family == "Greengenes_unannotated" && e.order == "Synechococcales" for e in values(gg2.taxonomy))
    unite = load_database(joinpath(TAXONOMY_DATA, "unite-10.0"))
    @test unite.taxonomy["SH0000001.10FU"].species == "Lactobacillus acidophilus"
    mktempdir() do directory
        # The marker file names the type of a directory with any name.
        custom = joinpath(directory, "custom")
        cp(joinpath(TAXONOMY_DATA, "unite-10.0"), custom)
        write(joinpath(custom, ".savont_db"), "unite-10.0\n")
        @test load_database(custom).taxonomy == unite.taxonomy
        @test_throws ArgumentError load_database(mkpath(joinpath(directory, "unknown")))
        # D10: a gzip-compressed SILVA taxonomy map.
        compressed = joinpath(directory, "silva-138.2")
        cp(joinpath(TAXONOMY_DATA, "silva-138.2"), compressed)
        taxmap = joinpath(compressed, "taxmap_slv_ssu_ref_nr_138.2.txt")
        write(taxmap * ".gz", transcode(GzipCompressor, read(taxmap))); rm(taxmap)
        @test Logging.with_logger(() -> load_database(compressed), Logging.NullLogger()).taxonomy == silva.taxonomy
        @test_throws ArgumentError load_database(compressed; upstream_quirks=true)
    end
end

@testset "classify, sintax and export" begin
    mktempdir() do directory
        quiet(f) = Logging.with_logger(f, Logging.NullLogger())
        single = joinpath(TAXONOMY_DATA, "asv_single")
        result = quiet(() -> savont_classify(single, joinpath(TAXONOMY_DATA, "emu-1"); output_dir=joinpath(directory, "c")))
        @test length(result.classifications) == 18
        @test issorted(result.classifications; by=c -> c.abundance, rev=true)
        @test count(c -> c.taxonomy === nothing, result.classifications) == 3   # no hit, missing taxonomy key, too short
        # EM abundances add up to the relative depth of the mapped ASVs (up to convergence).
        @test sum(result.abundances) ≈ sum(c.abundance for c in result.classifications if c.taxonomy !== nothing) atol = 1e-3
        # The tied ASV goes to the last of its equally abundant taxa, as upstream.
        tied = only(c for c in result.classifications if c.asv_header == "final_asv_8_depth_80")
        @test tied.taxonomy.species == "Roseburia intestinalis" && last(filter(m -> m.asv_header == tied.asv_header, result.mappings)).taxonomy == tied.taxonomy
        # EM moves the ASV equally close to two references to the one with more support.
        shared = only(c for c in result.classifications if c.asv_header == "final_asv_9_depth_70")
        @test shared.taxonomy.species == "Burkholderia mallei"
        hits = quiet(() -> savont_sintax(single, joinpath(TAXONOMY_DATA, "emu-1"); output_dir=joinpath(directory, "s")))
        @test length(hits) == 16 && all(h.bootstraps[2] > 0 for h in hits)   # all but the random and the too short ASV
        @test first(hits).names[2] == "Lactobacillus" && first(hits).bootstraps[2] == 1.0
        @test isfile(joinpath(directory, "s", "genus_abundance.tsv")) && !isfile(joinpath(directory, "s", "species_abundance.tsv"))
        # Fuzzy merging also absorbs ASVs trimmed at the 5' end, unlike upstream (D1).
        runs = joinpath.(TAXONOMY_DATA, ["export_a", "export_b", "export_c"])
        merged = quiet(() -> savont_export(runs, joinpath(directory, "e1")))
        upstream = quiet(() -> savont_export(runs, joinpath(directory, "e2"); upstream_quirks=true))
        unmerged = quiet(() -> savont_export(runs, joinpath(directory, "e3"); fuzzy=false))
        @test (length(merged.table), length(upstream.table), length(unmerged.table)) == (8, 10, 11)
        @test merged.sample_names == ["runA", "runB1", "runB2", "export_c"]
        @test sum(sum(last(v)) for v in values(merged.table)) == sum(sum(last(v)) for v in values(unmerged.table))
        @test_throws ArgumentError quiet(() -> savont_export(runs, joinpath(directory, "e4"); relabel=["x"]))
        @test (@test_logs (:warn, r"Duplicate sample names") match_mode=:any savont_export(runs, joinpath(directory, "e5");
                                                                                       relabel=["x", "x", "y", "z"])).sample_names[1] == "x"
        # D14: a rerun without ASVs leaves no result of the earlier run.
        empty_run = mkpath(joinpath(directory, "empty_run")); write(joinpath(empty_run, "final_asvs.fasta"), "")
        quiet(() -> savont_export([empty_run], joinpath(directory, "e1")))
        @test !isfile(joinpath(directory, "e1", "merged_taxon_counts.tsv"))
        @test length(readlines(joinpath(directory, "e1", "merged_feature_table.tsv"))) == 1
        @test isempty(quiet(() -> savont_sintax(empty_run, joinpath(TAXONOMY_DATA, "emu-1"); output_dir=joinpath(directory, "s"))))
        @test [length(readlines(joinpath(directory, "s", f))) for f in ("genus_abundance.tsv", "asv_mappings.tsv")] == [1, 1]
        # D15: a feature table without counts for the ASVs.
        zero_run = mkpath(joinpath(directory, "zero_run"))
        cp(joinpath(single, "final_asvs.fasta"), joinpath(zero_run, "final_asvs.fasta"))
        write(joinpath(zero_run, "feature-table.tsv"), "#OTU ID\tsampleA\nother_asv\t12\n")
        @test_logs (:warn, r"18 of 18 ASVs") match_mode=:any @test_throws ArgumentError savont_classify(zero_run, joinpath(TAXONOMY_DATA, "emu-1"))
        upstream_zero = quiet(() -> savont_classify(zero_run, joinpath(TAXONOMY_DATA, "emu-1"); upstream_quirks=true))
        @test all(c -> isnan(c.abundance), upstream_zero.classifications)
    end
end

@testset "failing tasks (SV6), pooled export without table rows (D16)" begin
    # Every task has finished when _tmap rethrows the first failure.
    finished = Threads.Atomic{Bool}(false)
    work(i) = i == 1 ? error("first task fails") : (sleep(0.5); finished[] = true)
    @test_throws TaskFailedException Savont._tmap(work, 1:2)
    @test finished[]
    @test Savont._tmap(x -> 2x, 1:3) == [2, 4, 6]

    mktempdir() do directory
        run = mkpath(joinpath(directory, "pooled"))
        sequence = String([b"ACGT"[(mm_hash64(UInt64(i)) & 3) + 1] for i in 1:300])
        write(joinpath(run, "final_asvs.fasta"), ">final_asv_0_depth_3-7\n$sequence\n")
        write(joinpath(run, "feature-table.tsv"), "#OTU ID\tS1\tS2\nfinal_asv_9_depth_1-1\t1\t1\n")
        counts(; kwargs...) = only(values(savont_export([run], joinpath(directory, "out"); fuzzy=false, kwargs...).table))[2]
        @test (@test_logs (:warn, r"1 ASVs .* not in its feature-table") counts()) == [3, 7]
        @test Logging.with_logger(() -> counts(upstream_quirks=true), Logging.NullLogger()) == [10, 0]
        write(joinpath(run, "final_asvs.fasta"), ">final_asv_0_depth_10\n$sequence\n")
        @test_throws ArgumentError Logging.with_logger(() -> counts(), Logging.NullLogger())
    end
end

@testset "command line" begin
    run_cli(args...) = (out = IOBuffer(); err = IOBuffer(); status = Savont.main(collect(String, args); output=out, err);
                        (status, String(take!(out)), String(take!(err))))
    @test run_cli("--help")[1] == 0 && occursin("classify", run_cli("--help")[2])
    @test occursin("--min-read-length", run_cli("asv", "--help")[2])
    @test occursin("savont 0.8.0", run_cli("--version")[2])
    @test run_cli()[1] == 2 && run_cli("cluster")[1] == 2 && run_cli("asv")[1] == 2
    @test run_cli("classify", "-i", "x")[1] == 2                        # --db required
    @test run_cli("asv", "r.fq", "--min-read-length", "ten")[1] == 2
    @test run_cli("asv", "r.fq", "--hifi=yes")[1] == 2
    @test run_cli("download", "--location", "x", "--dbs", "emu-2")[1] == 2
    command, values, positional = Savont.parse_command_line(["asv", "a.fq", "-m", "900", "--max-read-length=3000",
                                                             "b.fq", "-c7", "--rrna-operon"])
    @test (command, positional) == ("asv", ["a.fq", "b.fq"])
    config = Savont._asv_config(values)
    @test (config.min_read_length, config.max_read_length, config.sampling_rate) == (3500, 5000, 7)   # preset wins, as upstream
    @test Savont.parse_command_line(["export", "-i", "a", "b", "-o", "c", "--relabel", "x", "y"])[2]["input-dirs"] == ["a", "b"]
    mktempdir() do directory
        quiet(f) = Logging.with_logger(f, Logging.NullLogger())
        inputs = joinpath.(@__DIR__, "data", ["amplicons_sample1.fastq.gz", "amplicons_sample2.fastq.gz"])
        status, _, message = quiet(() -> run_cli("asv", inputs..., "--pooled-samples", "-o", joinpath(directory, "cli")))
        @test status == 0 && occursin("ASVs", message)
        quiet(() -> savont_asv(inputs, SavontConfig(pooled_samples=true); output_dir=joinpath(directory, "api")))
        for f in ("final_asvs.fasta", "feature-table.tsv", "final_assignments.tsv")
            @test read(joinpath(directory, "cli", f)) == read(joinpath(directory, "api", f))
        end
        @test run_cli("asv", inputs[1], "-c", "10", "-o", joinpath(directory, "even"))[1] == 1   # D11
        single = joinpath(TAXONOMY_DATA, "asv_single"); emu = joinpath(TAXONOMY_DATA, "emu-1")
        @test quiet(() -> run_cli("classify", "-i", single, "-d", emu, "-o", joinpath(directory, "c"), "--species-threshold", "99.8"))[1] == 0
        quiet(() -> savont_classify(single, emu; output_dir=joinpath(directory, "c_api"), species_threshold=99.8))
        @test read(joinpath(directory, "c", "species_abundance.tsv")) == read(joinpath(directory, "c_api", "species_abundance.tsv"))
        @test quiet(() -> run_cli("sintax", "-i", single, "-d", emu, "-o", joinpath(directory, "s"), "--n-iter", "20"))[1] == 0
        @test isfile(joinpath(directory, "s", "asv_mappings.tsv"))
        runs = joinpath.(TAXONOMY_DATA, ["export_a", "export_b"])
        @test quiet(() -> run_cli("export", "-i", runs..., "-o", joinpath(directory, "e"), "--relabel", "x", "y", "z"))[1] == 0
        @test startswith(readline(joinpath(directory, "e", "merged_feature_table.tsv")), "#OTU ID\tx\ty\tz")
        @test quiet(() -> run_cli("export", "-i", runs..., "-o", joinpath(directory, "e2"), "--relabel", "x"))[1] == 1
        @test quiet(() -> run_cli("classify", "-i", joinpath(directory, "missing"), "-d", emu))[1] == 1
    end
end

# Windows only: run `f` while `path` (a file or a directory) is open without
# FILE_SHARE_DELETE, as another program may hold it, so that it cannot be renamed.
function with_locked_path(f, path)
    handle = ccall((:CreateFileW, "kernel32"), stdcall, Ptr{Cvoid},
                   (Cwstring, UInt32, UInt32, Ptr{Cvoid}, UInt32, UInt32, Ptr{Cvoid}),
                   path, 0x80000000, 0x00000003, C_NULL, 3, 0x02000000, C_NULL)   # GENERIC_READ, share read/write, OPEN_EXISTING, BACKUP_SEMANTICS
    handle == Ptr{Cvoid}(-1 % UInt) && error("CreateFileW failed for $path")
    try
        return f()
    finally
        ccall((:CloseHandle, "kernel32"), stdcall, Int32, (Ptr{Cvoid},), handle)
    end
end

@testset "database download (local archives)" begin
    mktempdir() do directory
        source = mkpath(joinpath(directory, "source"))
        url(path) = "file:///" * replace(path, '\\' => '/')
        gzip(path) = (write(path * ".gz", transcode(GzipCompressor, read(path))); path * ".gz")
        function tgz(root, name)
            tarball = joinpath(source, name)
            open(tarball, "w") do io
                stream = GzipCompressorStream(io)
                Tar.create(root, stream)
                close(stream)
            end
            return tarball
        end
        emu_root = mkpath(joinpath(directory, "emu_tree", "emu_default"))
        foreach(f -> cp(joinpath(TAXONOMY_DATA, "emu-1", f), joinpath(emu_root, f)), ("species_taxid.fasta", "taxonomy.tsv"))
        emu_archive = tgz(dirname(emu_root), "emu.tar.gz")
        silva_files = [joinpath(TAXONOMY_DATA, "silva-138.2", "SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz"),
                       gzip(cp(joinpath(TAXONOMY_DATA, "silva-138.2", "taxmap_slv_ssu_ref_nr_138.2.txt"),
                               joinpath(source, "taxmap_slv_ssu_ref_nr_138.2.txt")))]
        # A raw UNITE FASTA: Species|Accession|SH|type|lineage headers, one malformed.
        unite_root = mkpath(joinpath(directory, "unite_tree"))
        open(joinpath(unite_root, "sh_general_release_dynamic_all.fasta"), "w") do io
            for (header, sequence) in read_fasta_records(joinpath(TAXONOMY_DATA, "unite-10.0", "unite_sequences.fasta"))
                lineage = get(Dict(eachline(joinpath(TAXONOMY_DATA, "unite-10.0", "unite_taxonomy.tsv")) .|>
                                   (l -> Pair(split(l, '\t')...))), header, nothing)
                lineage === nothing ? println(io, ">", header, "_", "x"^300) : println(io, ">Name|ACC1|", header, "|refs|", lineage)
                println(io, String(sequence))
            end
        end
        unite_archive = tgz(unite_root, "unite.tgz")
        urls = Dict("emu-1" => [url(emu_archive)], "silva-138.2" => url.(silva_files),
                    "greengenes2-2024.09" => [url(joinpath(TAXONOMY_DATA, "greengenes2-2024.09", "gg2_2024_09_toSpecies_trainset.fa.gz"))],
                    "unite-10.0" => [url(unite_archive)])
        location = joinpath(directory, "databases")
        quiet(f) = Logging.with_logger(f, Logging.NullLogger())
        directories = quiet(() -> download_database(location, "emu-1", "silva-138.2", "greengenes2-2024.09", "unite-10.0"; urls))
        @test basename.(directories) == ["emu-1", "silva-138.2", "greengenes2-2024.09", "unite-10.0"]
        for d in directories
            @test strip(read(joinpath(d, ".savont_db"), String)) == basename(d)
            renamed = joinpath(directory, "renamed_" * basename(d))
            cp(d, renamed)   # typed by the marker, not the name
            @test quiet(() -> load_database(renamed)).taxonomy == quiet(() -> load_database(joinpath(TAXONOMY_DATA, basename(d)))).taxonomy
        end
        @test sort(readdir(joinpath(location, "emu-1"))) == [".savont_db", "species_taxid.fasta", "taxonomy.tsv"]
        @test isfile(joinpath(location, "silva-138.2", "taxmap_slv_ssu_ref_nr_138.2.txt"))
        unite_names = first.(read_fasta_records(joinpath(location, "unite-10.0", "unite_sequences.fasta")))
        @test all(n -> startswith(n, "SH"), unite_names) && any(n -> length(n) == 255, unite_names)
        @test_throws ArgumentError download_database(location, "emu-2")

        # D20 (SV12): a download that fails halfway leaves the previous database whole.
        silva = joinpath(location, "silva-138.2")
        before = Dict(f => read(joinpath(silva, f)) for f in readdir(silva))
        new_fasta = joinpath(source, "new.fasta"); write(new_fasta, ">NEW1\nACGTACGTAC\n")
        broken = Dict("silva-138.2" => [url(gzip(new_fasta)), url(joinpath(source, "missing.txt.gz"))])
        @test_throws Exception quiet(() -> download_database(location, "silva-138.2"; urls=broken))
        @test Dict(f => read(joinpath(silva, f)) for f in readdir(silva)) == before
        @test sort(readdir(location)) == ["emu-1", "greengenes2-2024.09", "silva-138.2", "unite-10.0"]   # no staging left
        @test quiet(() -> load_database(silva)).taxonomy == quiet(() -> load_database(joinpath(TAXONOMY_DATA, "silva-138.2"))).taxonomy
        # D20 (SV13): a complete download with invalid content (an HTML page) does not replace it either.
        html = joinpath(source, "page.fasta"); write(html, "<!DOCTYPE html>\n<html><body>Not found</body></html>\n")
        invalid = Dict("silva-138.2" => [url(gzip(html)), url(silva_files[2])])
        @test_throws ArgumentError quiet(() -> download_database(location, "silva-138.2"; urls=invalid))
        @test Dict(f => read(joinpath(silva, f)) for f in readdir(silva)) == before
        @test sort(readdir(location)) == ["emu-1", "greengenes2-2024.09", "silva-138.2", "unite-10.0"]
        if Sys.iswindows()
            # SV16: the previous directory cannot be moved aside (open in another program):
            # nothing is copied or deleted, the download fails and the database is intact.
            @test_throws Base.IOError with_locked_path(() -> quiet(() -> download_database(location, "silva-138.2"; urls)), silva)
            @test Dict(f => read(joinpath(silva, f)) for f in readdir(silva)) == before
            @test sort(readdir(location)) == ["emu-1", "greengenes2-2024.09", "silva-138.2", "unite-10.0"]
        end
    end
end

@testset "input validation (D17, D18, D19, D22, SV10)" begin
    quiet(f) = Logging.with_logger(f, Logging.NullLogger())
    mktempdir() do directory
        sequence(seed) = String([b"ACGT"[(mm_hash64(UInt64(i + seed)) & 3) + 1] for i in 1:300])
        first_run = mkpath(joinpath(directory, "run1")); second_run = mkpath(joinpath(directory, "run2"))
        write(joinpath(first_run, "final_asvs.fasta"), ">final_asv_0_depth_7\n$(sequence(0))\n")
        write(joinpath(first_run, "feature-table.tsv"), "#OTU ID\tS1\nfinal_asv_0_depth_7\t7\n")
        write(joinpath(second_run, "feature-table.tsv"), "#OTU ID\tS2\nfinal_asv_0_depth_5\t5\n")
        output = joinpath(directory, "out")
        # D18 (SV8): a run without readable ASVs is an error, before any output is written.
        @test_throws ArgumentError quiet(() -> savont_export([first_run, second_run], output; fuzzy=false))
        @test !isdir(output) || isempty(readdir(output))
        counts = quiet(() -> savont_export([first_run, second_run], output; fuzzy=false, upstream_quirks=true)).table
        @test only(values(counts))[2] == [7, 0]   # upstream
        @test quiet(() -> Savont.main(["export", "-i", first_run, second_run, "-o", joinpath(directory, "cli")];
                                      output=IOBuffer(), err=IOBuffer())) == 1
        # D17 (SV9): counts that are not non-negative integers are errors, not zeros.
        write(joinpath(second_run, "final_asvs.fasta"), ">final_asv_0_depth_5\n$(sequence(1000))\n")
        for cell in ("5.0", "five", "-5", "")
            write(joinpath(second_run, "feature-table.tsv"), "#OTU ID\tS2\nfinal_asv_0_depth_5\t$cell\n")
            @test_throws ArgumentError quiet(() -> savont_export([first_run, second_run], output; fuzzy=false))
            @test_throws ArgumentError Savont._parse_feature_table(joinpath(second_run, "feature-table.tsv"))
            @test Savont._parse_feature_table(joinpath(second_run, "feature-table.tsv"); upstream_quirks=true)[2] ==
                  Dict("final_asv_0_depth_5" => [0])
        end
        write(joinpath(second_run, "feature-table.tsv"), "#OTU ID\tS2\tS3\nfinal_asv_0_depth_5\t5\n")   # missing cell
        @test_throws ArgumentError Savont._parse_feature_table(joinpath(second_run, "feature-table.tsv"))
        write(joinpath(second_run, "feature-table.tsv"), "#OTU ID\tS2\nfinal_asv_0_depth_5\t+5\n")
        @test sort([c[2] for c in values(quiet(() -> savont_export([first_run, second_run], output; fuzzy=false)).table)]) ==
              [[0, 5], [7, 0]]
        # D22 (SV18): a repeated OTU id is an error naming both lines, before any output.
        mktempdir() do run
            write(joinpath(run, "final_asvs.fasta"), ">final_asv_0_depth_12\n$(sequence(2000))\n")
            table = joinpath(run, "feature-table.tsv")
            for rows in (("7", "5"), ("5", "7"))
                write(table, "#OTU ID\tS1\nfinal_asv_0_depth_12\t$(rows[1])\nfinal_asv_0_depth_12\t$(rows[2])\n")
                error = try Savont._parse_feature_table(table); nothing catch err; err end
                @test error isa ArgumentError && occursin("lines 2 and 3", error.msg)
                @test Savont._parse_feature_table(table; upstream_quirks=true)[2] ==
                      Dict("final_asv_0_depth_12" => [parse(Int, rows[2])])   # upstream: the last row
                exported = joinpath(run, "out")
                @test_throws ArgumentError quiet(() -> savont_export([run], exported; fuzzy=false))
                @test !isdir(exported) || isempty(readdir(exported))
                @test quiet(() -> Savont.main(["export", "-i", run, "-o", joinpath(run, "cli")];
                                              output=IOBuffer(), err=IOBuffer())) == 1
                @test_throws ArgumentError Savont._read_feature_table(table, ["final_asv_0_depth_12"])
            end
            write(table, "#OTU ID\tS1\nfinal_asv_0_depth_12\t12\n")
            @test only(values(quiet(() -> savont_export([run], joinpath(run, "out"); fuzzy=false)).table))[2] == [12]
        end

        # D19 (SV11): preparing UNITE never truncates its source or earlier outputs.
        unite = mkpath(joinpath(directory, "unite"))
        aliased = joinpath(unite, "unite_sequences.fasta")
        write(aliased, ">Name|ACC1|SH1.10FU|refs|k__Fungi;s__X\nACGT\n")
        preprocess_unite(aliased, unite)
        @test read(aliased, String) == ">SH1.10FU\nACGT\n"
        @test readlines(joinpath(unite, "unite_taxonomy.tsv")) == ["sh_id\ttaxonomy", "SH1.10FU\tk__Fungi;s__X"]
        previous = Dict(f => read(joinpath(unite, f)) for f in readdir(unite))
        @test_throws ArgumentError preprocess_unite(joinpath(unite, "missing.fasta"), unite)
        @test Dict(f => read(joinpath(unite, f)) for f in readdir(unite)) == previous
        # SV15: an output path that is a directory is refused, not deleted.
        blocked = mkpath(joinpath(directory, "blocked"))
        mkpath(joinpath(blocked, "unite_sequences.fasta")); write(joinpath(blocked, "unite_sequences.fasta", "unrelated.txt"), "keep")
        @test_throws ArgumentError preprocess_unite(joinpath(unite, "unite_sequences.fasta"), blocked)
        @test read(joinpath(blocked, "unite_sequences.fasta", "unrelated.txt"), String) == "keep"
        @test sort(readdir(blocked)) == ["unite_sequences.fasta"]   # no temporary left
        # SV14: if publishing the second file fails, both previous files are restored.
        new_fasta = joinpath(unite, "new_fasta.tmp"); write(new_fasta, ">SH2\nTTTT\n")
        targets = joinpath.(unite, ["unite_sequences.fasta", "unite_taxonomy.tsv"])
        @test_throws Exception Savont._publish_files([new_fasta, joinpath(unite, "missing.tmp")], targets)
        @test Dict(f => read(joinpath(unite, f)) for f in readdir(unite) if f != "new_fasta.tmp") == previous
        if Sys.iswindows()
            # SV14 with a real lock on the taxonomy: the pair is kept and no backup copy is left.
            new_taxonomy = joinpath(unite, "new_taxonomy.tmp"); write(new_taxonomy, "sh_id\ttaxonomy\nSH2\tk__Fungi\n")
            @test_throws Base.IOError with_locked_path(() -> Savont._publish_files([new_fasta, new_taxonomy], targets), targets[2])
            @test Dict(f => read(joinpath(unite, f)) for f in readdir(unite) if !endswith(f, ".tmp")) == previous
            @test !any(f -> occursin(".previous-", f), readdir(unite))
        end
    end
    # SV10: qualities outside 0..=255 are refused, as clap does for a u8.
    for quality in ("256", "300")
        options = Savont.parse_command_line(["asv", "r.fq", "--minimum-base-quality=$quality"])[2]
        @test_throws Savont.CLIError Savont._asv_config(options)
    end
    options = Savont.parse_command_line(["asv", "r.fq", "--minimum-base-quality", "255"])[2]
    @test Savont._asv_config(options).minimum_base_quality == 0xff
end

include("rust_oracle.jl")
