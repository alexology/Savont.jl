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

function main()
    directory = joinpath(@__DIR__, "..", "test", "data")
    mkpath(directory)
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
end

main()
