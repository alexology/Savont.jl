# Generate the deterministic synthetic inputs in test/data.
#
#   julia --project=. tools/make_test_data.jl
#
# A SplitMix64 generator keeps the output independent of the Julia version.

using CodecZlib

mutable struct SplitMix64
    state::UInt64
end
function next!(rng::SplitMix64)
    rng.state += 0x9e3779b97f4a7c15
    z = rng.state
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    return z ⊻ (z >> 31)
end
uniform(rng) = (next!(rng) >> 11) * 2.0^-53
randint(rng, n) = Int(next!(rng) % UInt64(n)) + 1

const DNA = b"ACGT"
complement(base) = base == UInt8('A') ? UInt8('T') : base == UInt8('T') ? UInt8('A') :
                   base == UInt8('C') ? UInt8('G') : base == UInt8('G') ? UInt8('C') : UInt8('N')
revcomp(sequence) = reverse!(map(complement, sequence))

function mutate(rng, sequence, positions)
    result = copy(sequence)
    for position in positions
        result[position] = DNA[mod(findfirst(==(result[position]), DNA) + randint(rng, 3) - 1, 4) + 1]
    end
    return result
end

# Sequencing errors with matching low qualities; returns (sequence, qualities).
function sequence_read(rng, template; substitution=0.004, deletion=0.001, insertion=0.001)
    sequence = UInt8[]; qualities = UInt8[]
    for base in template
        r = uniform(rng)
        if r < deletion
            continue
        elseif r < deletion + insertion
            push!(sequence, DNA[randint(rng, 4)]); push!(qualities, UInt8(33 + 8 + randint(rng, 6)))
            push!(sequence, base); push!(qualities, UInt8(33 + 30 + randint(rng, 10)))
        elseif r < deletion + insertion + substitution
            push!(sequence, mutate(rng, [base], 1:1)[1]); push!(qualities, UInt8(33 + 5 + randint(rng, 15)))
        else
            push!(sequence, base); push!(qualities, UInt8(33 + 28 + randint(rng, 14)))
        end
    end
    return sequence, qualities
end

function write_fastq(io, name, sequence, qualities)
    write(io, '@', name, '\n', sequence, "\n+\n", qualities, '\n')
end

function amplicon_reads(path, rng, haplotypes, abundances, sample)
    open(path, "w") do raw
        io = endswith(path, ".gz") ? GzipCompressorStream(raw) : raw
        number = 0
        for (haplotype, abundance) in zip(haplotypes, abundances), _ in 1:abundance
            number += 1
            # Reads start and end up to 5 bases inside the amplicon.
            first_base = randint(rng, 6); last_base = length(haplotype) - randint(rng, 6) + 1
            sequence, qualities = sequence_read(rng, haplotype[first_base:last_base])
            name = "s$(sample)_read$(number)"
            r = uniform(rng)
            if r < 0.45
                sequence = revcomp(sequence); reverse!(qualities)
            elseif r < 0.55
                # Cutadapt-style mark: stored reverse-complemented, flipped back on input.
                sequence = revcomp(sequence); reverse!(qualities)
                name *= " rc"
            end
            uniform(rng) < 0.05 && (sequence[randint(rng, length(sequence))] = UInt8('N'))
            write_fastq(io, name, sequence, qualities)
        end
        # Edge cases: reads below the length filter, low estimated identity and
        # constant qualities (identity not estimated).
        for i in 1:2
            sequence, qualities = sequence_read(rng, haplotypes[1][1:600])
            write_fastq(io, "s$(sample)_short$i", sequence, qualities)
            sequence, _ = sequence_read(rng, haplotypes[1])
            write_fastq(io, "s$(sample)_lowq$i", sequence, [UInt8(33 + 9 + randint(rng, 6)) for _ in sequence])
        end
        sequence, _ = sequence_read(rng, haplotypes[1])
        write_fastq(io, "s$(sample)_constq", sequence, fill(UInt8('I'), length(sequence)))
        io === raw || close(io)
    end
end

# Short random records over an extended alphabet for primitive comparisons.
# The FASTA file repeats the sequences without qualities.
function primitive_records(path, rng, count)
    alphabet = b"ACGTACGTACGTacgtNRU"
    fasta = open(replace(path, ".fastq" => ".fasta"), "w")
    open(path, "w") do io
        for i in 1:count
            n = randint(rng, 80)
            sequence = [alphabet[randint(rng, length(alphabet))] for _ in 1:n]
            constant = uniform(rng) < 0.15
            qualities = constant ? fill(UInt8(33 + randint(rng, 40)), n) : [UInt8(33 + randint(rng, 40)) for _ in 1:n]
            write_fastq(io, "p$i", sequence, qualities)
            write(fasta, ">p", string(i), '
', sequence, '
')
        end
    end
    close(fasta)
end

# Closely related haplotypes (1, 2 and 7 SNPs from the first) with a higher error
# rate, so that greedy SNPmer clusters need merging and reassignment.
function noisy_reads(path, rng, substitution)
    base = [DNA[randint(rng, 4)] for _ in 1:1200]
    haplotypes = [base, mutate(rng, base, [300]), mutate(rng, base, [500, 900]), mutate(rng, base, 200:150:1100)]
    open(path, "w") do raw
        io = GzipCompressorStream(raw)
        number = 0
        for (haplotype, abundance) in zip(haplotypes, (220, 60, 25, 14)), _ in 1:abundance
            number += 1
            sequence, qualities = sequence_read(rng, haplotype; substitution, deletion=0.003, insertion=0.003)
            uniform(rng) < 0.5 && (sequence = revcomp(sequence); reverse!(qualities))
            write_fastq(io, "n$number", sequence, qualities)
        end
        close(io)
    end
end

# Reads of one haplotype for partial order alignment. Every `every`-th read moves
# `shift` bases from near the start to near the end (length unchanged), so the
# banded alignment misses the optimum and differs from the exact one.
function poa_reads(path, rng, reads, len, substitution; shift=0, every=0)
    haplotype = [DNA[randint(rng, 4)] for _ in 1:len]
    open(path, "w") do raw
        io = GzipCompressorStream(raw)
        for r in 1:reads
            template = haplotype[randint(rng, 8):(len - randint(rng, 8) + 1)]
            if every > 0 && r % every == 0
                a = 150 + randint(rng, 50)
                template = vcat(haplotype[1:a], haplotype[(a + shift + 1):(end - 200)],
                                [DNA[randint(rng, 4)] for _ in 1:shift], haplotype[(end - 199):end])
            end
            sequence, qualities = sequence_read(rng, template; substitution,
                                                deletion=substitution / 2, insertion=substitution / 2)
            write_fastq(io, "r$r", sequence, qualities)
        end
        close(io)
    end
end

# Clusters for stage 4b, run with --low-polymorphism so that stage 3 does not
# split haplotypes: a 60/40 mix of two haplotypes differing at three central sites
# and at one site within 100 bases of each end (low-confidence positions, a
# low-quality consensus and extended end masks); a clean amplicon whose reads
# often start or end inside it (coverage trimming); a 96-base amplicon with a
# mixed site (consensus shorter than 100 positions); homopolymer runs throughout.
function mixed_reads(path, rng)
    with_runs(n) = vcat([vcat([DNA[randint(rng, 4)] for _ in 1:20], fill(DNA[randint(rng, 4)], 2 + randint(rng, 4))) for _ in 1:cld(n, 24)]...)[1:n]
    mixed = with_runs(1200)
    other = mutate(rng, mixed, [70, 400, 650, 900, 1150])
    clean = with_runs(1000)
    short = with_runs(96)
    short_other = mutate(rng, short, [48])
    groups = [(mixed, 36, 0), (other, 24, 0), (clean, 60, 120), (short, 14, 0), (short_other, 10, 0)]
    open(path, "w") do raw
        io = GzipCompressorStream(raw)
        number = 0
        for (haplotype, abundance, wander) in groups, _ in 1:abundance
            number += 1
            first_base = randint(rng, 4); last_base = length(haplotype) - randint(rng, 4) + 1
            if wander > 0 && uniform(rng) < 0.4
                first_base += randint(rng, wander); last_base -= randint(rng, wander)
            end
            sequence, qualities = sequence_read(rng, haplotype[first_base:last_base]; substitution=0.01, deletion=0.004, insertion=0.004)
            uniform(rng) < 0.5 && (sequence = revcomp(sequence); reverse!(qualities))
            write_fastq(io, "m$number", sequence, qualities)
        end
        close(io)
    end
end

# Chimeras for stage 6: (a) a chimera of two parents differing at 8% of sites;
# (b) a chimera of two parents differing at 1.5%, too shallow (less than 10 times
# the chimera) for upstream's parent-similarity exception; (c) a shallow variant of
# a deep consensus with 5 substitutions 40-100 bases from its end (single-parent
# rule). Run with the default settings.
function chimera_reads(path, rng)
    a = [DNA[randint(rng, 4)] for _ in 1:1300]
    b = mutate(rng, a, 7:13:1300)
    c = [DNA[randint(rng, 4)] for _ in 1:1300]
    d = mutate(rng, c, 30:65:1300)
    e = [DNA[randint(rng, 4)] for _ in 1:1300]
    e_variant = mutate(rng, e, 1200:12:1260)
    groups = [(a, 150), (b, 100), (vcat(a[1:650], b[651:end]), 15), (c, 140), (d, 130), (vcat(c[1:650], d[651:end]), 16),
              (e, 240), (e_variant, 12)]
    open(path, "w") do raw
        io = GzipCompressorStream(raw)
        number = 0
        for (haplotype, abundance) in groups, _ in 1:abundance
            number += 1
            first_base = randint(rng, 4); last_base = length(haplotype) - randint(rng, 4) + 1
            sequence, qualities = sequence_read(rng, haplotype[first_base:last_base])
            uniform(rng) < 0.5 && (sequence = revcomp(sequence); reverse!(qualities))
            write_fastq(io, "c$number", sequence, qualities)
        end
        close(io)
    end
end

# Taxonomy inputs for classify, sintax and export (test/data/taxonomy): ASV
# directories as written by `savont asv` and one small reference database per
# upstream format, all holding the same reference sequences.
const LINEAGES = [  # superkingdom, phylum, class, order, family, genus, species, clade
    ("Bacteria", "Bacillota", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus", "Lactobacillus acidophilus", "Terrabacteria group"),
    ("Bacteria", "Bacillota", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus", "Lactobacillus gasseri", "Terrabacteria group"),
    ("Bacteria", "Bacillota", "Bacilli", "Bacillales", "Bacillaceae", "Bacillus", "Bacillus subtilis", "Terrabacteria group"),
    ("Bacteria", "Pseudomonadota", "Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae", "Escherichia", "Escherichia coli", ""),
    ("Bacteria", "Pseudomonadota", "Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae", "Salmonella", "Salmonella enterica", ""),
    ("Bacteria", "Pseudomonadota", "Gammaproteobacteria", "Pseudomonadales", "Pseudomonadaceae", "Pseudomonas", "Pseudomonas aeruginosa", ""),
    ("Bacteria", "Actinomycetota", "Actinomycetes", "Mycobacteriales", "Mycobacteriaceae", "Mycobacterium", "Mycobacterium tuberculosis", "Terrabacteria group"),
    ("Archaea", "Methanobacteriota", "Methanobacteria", "Methanobacteriales", "Methanobacteriaceae", "Methanobrevibacter", "Methanobrevibacter smithii", ""),
    ("Bacteria", "Bacteroidota", "Bacteroidia", "Bacteroidales", "Bacteroidaceae", "Bacteroides", "Bacteroides fragilis", "FCB group"),
    ("Bacteria", "Bacillota", "Clostridia", "Eubacteriales", "Lachnospiraceae", "Roseburia", "Roseburia intestinalis", "Terrabacteria group"),
    ("Bacteria", "Pseudomonadota", "Betaproteobacteria", "Burkholderiales", "Burkholderiaceae", "Burkholderia", "Burkholderia cepacia", ""),
    ("Bacteria", "Pseudomonadota", "Betaproteobacteria", "Burkholderiales", "Burkholderiaceae", "Burkholderia", "Burkholderia mallei", ""),
    ("Bacteria", "Verrucomicrobiota", "Verrucomicrobiia", "Verrucomicrobiales", "Akkermansiaceae", "Akkermansia", "Akkermansia muciniphila", "PVC group"),
    ("Bacteria", "Bacillota", "Bacilli", "Lactobacillales", "Streptococcaceae", "Streptococcus", "Streptococcus mutans", "Terrabacteria group"),
    ("Bacteria", "Bacillota", "Bacilli", "Lactobacillales", "Enterococcaceae", "Enterococcus", "Enterococcus faecalis", "Terrabacteria group"),
    ("Bacteria", "Cyanobacteriota", "Cyanophyceae", "Synechococcales", "", "", "", "Terrabacteria group"),
]

random_dna(rng, n) = [DNA[randint(rng, 4)] for _ in 1:n]
# `count` substitutions at distinct random positions.
function substitute(rng, sequence, count)
    positions = Int[]
    while length(positions) < count
        p = randint(rng, length(sequence))
        p in positions || push!(positions, p)
    end
    return mutate(rng, sequence, positions)
end

# About `fraction` substitutions, confined to the second half of each 50-base
# block so that minimap2 still finds seeds at low identity.
function substitute_blocks(rng, sequence, fraction)
    result = copy(sequence)
    for start in 1:50:(length(sequence) - 49)
        positions = Int[]
        while length(positions) < round(Int, 50fraction)
            p = start + 24 + randint(rng, 25)
            p in positions || push!(positions, p)
        end
        result = mutate(rng, result, positions)
    end
    return result
end

wrap(sequence, width) = join((String(sequence[i:min(i + width - 1, end)]) for i in 1:width:length(sequence)), '\n')

function write_gzip(path, text)
    open(path, "w") do raw
        io = GzipCompressorStream(raw)
        write(io, text)
        close(io)
    end
end

# ASVs (sequence, depth) and references (sequence, lineage index or 0 for a key
# absent from the taxonomy), chosen to reach every identity band of
# `from_taxonomy_entry`, ties between identical references, EM redistribution,
# reverse strands, ambiguous and lower-case reference bases and unmapped ASVs.
function taxonomy_sets(rng)
    asvs = Tuple{Vector{UInt8},Int}[]
    references = Tuple{Vector{UInt8},Int}[]
    for (substitutions, lineage, depth) in ((0, 1, 500), (5, 2, 320), (30, 3, 210), (100, 4, 150),
                                           (160, 5, 120), (220, 6, 100), (235, 7, 95))
        reference = random_dna(rng, 1000 + randint(rng, 100))
        push!(references, (reference, lineage))
        mutated = substitutions < 150 ? substitute(rng, reference, substitutions) :
                  substitute_blocks(rng, reference, substitutions / 1000)
        push!(asvs, (mutated, depth))
    end
    references[2][1][400] = UInt8('N')
    references[3][1][1:120] .|= 0x20
    push!(asvs, (random_dna(rng, 1050), 90))                       # no hit
    tie = random_dna(rng, 1020)                                    # identical references
    push!(references, (tie, 9)); push!(references, (copy(tie), 10))
    push!(asvs, (substitute(rng, tie, 3), 80))
    shared = random_dna(rng, 1040)                                 # EM redistribution
    other = copy(shared); other[300] = mutate(rng, other[300:300], 1:1)[1]; other[700] = mutate(rng, other[700:700], 1:1)[1]
    push!(references, (shared, 11)); push!(references, (other, 12))
    between = copy(shared); between[700] = other[700]
    push!(asvs, (between, 70)); push!(asvs, (copy(other), 60))
    reverse_reference = random_dna(rng, 1010)
    push!(references, (reverse_reference, 13))
    push!(asvs, (revcomp(substitute(rng, reverse_reference, 2)), 50))
    missing_key = random_dna(rng, 1000)
    push!(references, (missing_key, 0))
    push!(asvs, (copy(missing_key), 40))
    push!(asvs, (random_dna(rng, 10), 30))                         # too short
    for lineage in (14, 15)                                        # equal abundances
        reference = random_dna(rng, 1000)
        push!(references, (reference, lineage))
        push!(asvs, (substitute(rng, reference, 1), 25))
    end
    push!(asvs, (substitute(rng, references[1][1], 1), 20))       # same species as ASV 0
    unannotated = random_dna(rng, 1000)
    push!(references, (unannotated, 16))
    push!(asvs, (copy(unannotated), 15))
    for i in 1:15                                                  # background
        push!(references, (random_dna(rng, 1000), mod1(i, 15)))
    end
    return asvs, references
end

asv_header(i, depth) = "final_asv_$(i - 1)_depth_$(depth) debug_id:$(i - 1) chimera_score:0 unambiguous_read_assignments:$(depth) ambig_read_assignments:0 num_align_leq_10_mismatches:$(depth)"

function write_asv_directory(directory, asvs; samples=nothing, feature_table=true)
    mkpath(directory)
    depth_field(depth, i) = samples === nothing ? string(depth) : join(samples[i], '-')
    open(joinpath(directory, "final_asvs.fasta"), "w") do io
        for (i, (sequence, depth)) in enumerate(asvs)
            write(io, '>', asv_header(i, depth_field(depth, i)), '\n', sequence, '\n')
        end
    end
    feature_table || return
    open(joinpath(directory, "feature-table.tsv"), "w") do io
        names = samples === nothing ? ["sampleA"] : ["S$(k)" for k in eachindex(first(samples))]
        println(io, "#OTU ID\t", join(names, '\t'))
        for (i, (_, depth)) in enumerate(asvs)
            counts = samples === nothing ? [depth] : samples[i]
            println(io, "final_asv_$(i - 1)_depth_$(depth_field(depth, i))\t", join(counts, '\t'))
        end
    end
end

function write_databases(directory, references)
    taxid(i) = string(1000 + i)
    # EMU: species_taxid.fasta and taxonomy.tsv (12 columns, one malformed line).
    emu = joinpath(directory, "emu-1"); mkpath(emu)
    open(joinpath(emu, "species_taxid.fasta"), "w") do io
        for (r, (sequence, lineage)) in enumerate(references)
            key = lineage == 0 ? "99999" : taxid(lineage)
            write(io, '>', key, ":emu_db:", string(r), " [reference $r]\n", wrap(sequence, 80), '\n')
        end
    end
    open(joinpath(emu, "taxonomy.tsv"), "w") do io
        println(io, "tax_id\tspecies\tgenus\tfamily\torder\tclass\tphylum\tclade\tsuperkingdom\tsubspecies\tspecies subgroup\tspecies group")
        for (l, (sk, p, c, o, f, g, s, clade)) in enumerate(LINEAGES)
            println(io, join([taxid(l), s, g, f, o, c, p, clade, sk, "", "", ""], '\t'))
        end
        println(io, "12345\tTruncated line")
    end
    # SILVA: RNA sequences in a gzip-compressed FASTA, taxmap with variable depth.
    silva = joinpath(directory, "silva-138.2"); mkpath(silva)
    accession(r) = "AB" * lpad(string(r), 6, '0')
    fasta = IOBuffer()
    for (r, (sequence, lineage)) in enumerate(references)
        rna = replace(String(copy(sequence)), 'T' => 'U', 't' => 'u')
        path = lineage == 0 ? "Bacteria;Unknown" : join(LINEAGES[lineage][1:7], ';')
        write(fasta, '>', accession(r), ".1.", string(length(sequence)), ' ', path, '\n', wrap(Vector{UInt8}(rna), 60), '\n')
    end
    write_gzip(joinpath(silva, "SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz"), take!(fasta))
    open(joinpath(silva, "taxmap_slv_ssu_ref_nr_138.2.txt"), "w") do io
        println(io, "primaryAccession\tstart\tstop\tpath\torganism_name\ttaxid")
        for (r, (sequence, lineage)) in enumerate(references)
            lineage == 0 && continue
            sk, p, c, o, f, g, s, _ = LINEAGES[lineage]
            # Lineage 16 stops at the order, as some SILVA paths do.
            path = lineage == 16 ? "$sk;$p;$c;$o;" : "$sk;$p;$c;$o;$f;$g;"
            println(io, join([accession(r), "1", string(length(sequence)), path, s, string(40000 + lineage)], '\t'))
        end
        println(io, "AB999999\t1\t100")
    end
    # GreenGenes2: the header is the taxonomy; lineage 1 has a space (minimap2
    # truncates the target name there) and lineage 16 empty ranks.
    gg2 = joinpath(directory, "greengenes2-2024.09"); mkpath(gg2)
    fasta = IOBuffer()
    for (sequence, lineage) in references
        lineage == 0 && continue
        sk, p, c, o, f, g, s, _ = LINEAGES[lineage]
        epithet = lineage == 1 ? "acidophilus group" : isempty(s) ? "" : last(split(s))
        write(fasta, ">d__$sk;p__$p;c__$c;o__$o;f__$f;g__$g;s__$epithet;\n", wrap(sequence, 70), '\n')
    end
    write_gzip(joinpath(gg2, "gg2_2024_09_toSpecies_trainset.fa.gz"), take!(fasta))
    # UNITE, as written by `savont download`: SH accessions and a taxonomy table.
    unite = joinpath(directory, "unite-10.0"); mkpath(unite)
    sh(r) = "SH" * lpad(string(r), 7, '0') * ".10FU"
    open(joinpath(unite, "unite_sequences.fasta"), "w") do io
        for (r, (sequence, _)) in enumerate(references)
            write(io, '>', sh(r), '\n', sequence, '\n')
        end
    end
    open(joinpath(unite, "unite_taxonomy.tsv"), "w") do io
        println(io, "sh_id\ttaxonomy")
        for (r, (_, lineage)) in enumerate(references)
            lineage == 0 && continue
            sk, p, c, o, f, g, s, _ = LINEAGES[lineage]
            println(io, sh(r), "\tk__$sk;p__$p;c__$c;o__$o;f__$f;g__$g;s__", replace(s, ' ' => '_'))
        end
    end
end

# Directories for `savont export`: a single-sample run with classify mappings, a
# pooled run with sintax mappings holding a reverse complement, 3'- and 5'-trimmed
# copies (fuzzy merging, D1) and a copy trimmed beyond 10 bp, and a run with only
# final_asvs.fasta.
function write_export_directories(directory, rng)
    base = [random_dna(rng, 1000 + 20i) for i in 1:6]
    a = [(base[1], 300), (base[2], 200), (base[3], 150), (base[4], 90)]
    write_asv_directory(joinpath(directory, "export_a"), a)
    mv(joinpath(directory, "export_a", "feature-table.tsv"), joinpath(directory, "export_a", "ft.tmp"))
    open(joinpath(directory, "export_a", "feature-table.tsv"), "w") do io
        println(io, "# Constructed from savont")
        write(io, replace(read(joinpath(directory, "export_a", "ft.tmp"), String), "sampleA" => "runA"))
    end
    rm(joinpath(directory, "export_a", "ft.tmp"))
    columns = "asv_header\tdepth\talignment_identity\tnumber_mismatches\ttax_id\tspecies\tgenus\tfamily\torder\tclass\tphylum\tclade\tsuperkingdom\treference"
    open(joinpath(directory, "export_a", "asv_mappings.tsv"), "w") do io
        println(io, columns)
        for (i, lineage, identity) in ((1, 1, "100.00"), (1, 2, "100.00"), (2, 4, "99.20"), (3, 9, "97.10"))
            sk, p, c, o, f, g, s, clade = LINEAGES[lineage]
            println(io, join(["final_asv_$(i - 1)_depth_$(a[i][2])", a[i][2], identity, "0", string(1000 + lineage),
                              s, g, f, o, c, p, clade, sk, "$(1000 + lineage):emu_db:$i"], '\t'))
        end
    end
    b_samples = [[60, 40], [7, 3], [12, 0], [5, 4], [0, 33], [21, 21]]
    b = [(revcomp(base[1]), 100), (base[2][1:(end - 5)], 10), (base[3][5:end], 12), (base[4][13:end], 9),
         (base[5], 33), (base[6], 42)]
    write_asv_directory(joinpath(directory, "export_b"), b; samples=b_samples)
    open(joinpath(directory, "export_b", "feature-table.tsv"), "w") do io
        println(io, "#OTU ID\trunB1\trunB2")
        for (i, counts) in enumerate(b_samples)
            println(io, "final_asv_$(i - 1)_depth_$(join(counts, '-'))\t", join(counts, '\t'))
        end
    end
    open(joinpath(directory, "export_b", "asv_mappings.tsv"), "w") do io
        println(io, "asv_header\tdepth\tspecies_bootstrap\tgenus_bootstrap\tfamily_bootstrap\torder_bootstrap\tclass_bootstrap\tphylum_bootstrap\tsuperkingdom_bootstrap\tspecies\tgenus\tfamily\torder\tclass\tphylum\tsuperkingdom")
        for (i, lineage) in ((1, 3), (3, 5), (5, 6), (6, 7))
            sk, p, c, o, f, g, _, _ = LINEAGES[lineage]
            println(io, join(["final_asv_$(i - 1)_depth_$(join(b_samples[i], '-'))", sum(b_samples[i]),
                              "0.500", "0.900", "0.950", "1.000", "1.000", "1.000", "1.000",
                              "UNCLASSIFIED", g, f, o, c, p, sk], '\t'))
        end
    end
    c = [(base[5], 17), (base[6][3:end], 8), (random_dna(rng, 900), 5)]
    write_asv_directory(joinpath(directory, "export_c"), c; feature_table=false)
    text = read(joinpath(directory, "export_c", "final_asvs.fasta"), String)
    write(joinpath(directory, "export_c", "final_asvs.fasta"), replace(text, "final_asv_2_depth_5 " => "final_asv_2_depth_x "))
end

function taxonomy_data(directory)
    rng = SplitMix64(0x7a8017)
    directory = joinpath(directory, "taxonomy")
    rm(directory; force=true, recursive=true)
    asvs, references = taxonomy_sets(rng)
    write_asv_directory(joinpath(directory, "asv_single"), asvs)
    samples = [[d ÷ 2, d - d ÷ 2 - d ÷ 5, d ÷ 5] for (_, d) in asvs]
    samples[3] = [0, 210, 0]; samples[9] = [80, 0, 0]
    write_asv_directory(joinpath(directory, "asv_pooled"), [(s, sum(n)) for ((s, _), n) in zip(asvs, samples)]; samples)
    write_asv_directory(joinpath(directory, "asv_headers"), asvs[1:2:end]; feature_table=false)
    write_databases(directory, references)
    write_export_directories(directory, rng)
end

function main()
    directory = joinpath(@__DIR__, "..", "test", "data")
    mkpath(directory)
    taxonomy_data(directory)
    "taxonomy" in ARGS && return
    rng = SplitMix64(0x5a7017)
    base = [DNA[randint(rng, 4)] for _ in 1:1300]
    variant_b = mutate(rng, base, 100:200:1300)
    variant_c = mutate(rng, variant_b, 150:180:1300)
    # Two unrelated amplicons; the last one is too rare to form a cluster.
    other = [DNA[randint(rng, 4)] for _ in 1:1250]
    rare = [DNA[randint(rng, 4)] for _ in 1:1400]
    haplotypes = [base, variant_b, variant_c, other, rare]
    amplicon_reads(joinpath(directory, "amplicons_sample1.fastq.gz"), rng, haplotypes, (70, 35, 15, 25, 4), 1)
    amplicon_reads(joinpath(directory, "amplicons_sample2.fastq.gz"), rng, haplotypes, (40, 40, 10, 20, 0), 2)
    primitive_records(joinpath(directory, "primitives.fastq"), rng, 60)
    noisy_reads(joinpath(directory, "noisy_1pct.fastq.gz"), SplitMix64(0x0bad5eed), 0.01)
    noisy_reads(joinpath(directory, "noisy_2pct.fastq.gz"), SplitMix64(0x0bad5eed), 0.02)
    poa_reads(joinpath(directory, "poa_plain.fastq.gz"), SplitMix64(0x9a0001), 76, 1300, 0.01)
    poa_reads(joinpath(directory, "poa_shifted.fastq.gz"), SplitMix64(0x9a0002), 40, 1300, 0.01; shift=250, every=3)
    poa_reads(joinpath(directory, "poa_shifted_short.fastq.gz"), SplitMix64(0x9a0003), 60, 900, 0.01; shift=120, every=2)
    mixed_reads(joinpath(directory, "mixed.fastq.gz"), SplitMix64(0x4b1e7d))
    chimera_reads(joinpath(directory, "chimeras.fastq.gz"), SplitMix64(0xc41e2a))
end

main()
