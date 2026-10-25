# The `savont asv` pipeline (upstream `main.rs::run_cluster`) and its output files.

const ASV_FILE = "final_asvs.fasta"

_depth_field(c::ConsensusSequence) = isempty(c.per_sample_depths) ? string(c.depth) : join(c.per_sample_depths, '-')

# Decompressed sequence without leading and trailing N (upstream
# `write_consensus_fasta`).
function _output_sequence(c::ConsensusSequence)
    sequence = homopolymer_decompress(c.sequence, c.hp_lengths)
    first = something(findfirst(!=(UInt8('N')), sequence), 1)
    last = something(findlast(!=(UInt8('N')), sequence), length(sequence))
    return sequence[first:last]
end

"""
    write_asv_fasta(path, consensuses, prefix)

Upstream `write_consensus_fasta`: one record per consensus, named
`<prefix>_asv_<i>_depth_<depth>` (0-based `i`; per-sample depths joined by `-`)
with the 0-based cluster index (`debug_id`), chimera score and assignment counts
in the description, and the decompressed sequence without masked ends.
"""
function write_asv_fasta(path::AbstractString, consensuses::AbstractVector{ConsensusSequence}, prefix::AbstractString)
    open(path, "w") do io
        for (i, c) in enumerate(consensuses)
            println(io, ">", prefix, "_asv_", i - 1, "_depth_", _depth_field(c), " debug_id:", c.id - 1,
                    " chimera_score:", something(c.chimera_score, 0),
                    " unambiguous_read_assignments:", something(c.unambiguous_count, 0),
                    " ambig_read_assignments:", something(c.ambiguous_count, 0),
                    " num_align_leq_10_mismatches:", something(c.close_count, 0))
            write(io, _output_sequence(c), '\n')
        end
    end
    return path
end

"""
    write_feature_table(path, consensuses, sample_names)

QIIME2-style table: `#OTU ID` and the sample names, then one row per ASV
(`final_asv_<i>_depth_<depth>` and its depth, or its per-sample depths).
"""
function write_feature_table(path::AbstractString, consensuses::AbstractVector{ConsensusSequence},
                             sample_names::AbstractVector{<:AbstractString})
    open(path, "w") do io
        println(io, "#OTU ID\t", join(sample_names, '\t'))
        for (i, c) in enumerate(consensuses)
            if isempty(c.per_sample_depths)
                println(io, "final_asv_", i - 1, "_depth_", c.depth, '\t', c.depth)
            else
                println(io, "final_asv_", i - 1, "_depth_", _depth_field(c), '\t', join(c.per_sample_depths, '\t'))
            end
        end
    end
    return path
end

"""
    write_final_assignments(path, reads, asvs, mappings; upstream_quirks=false)

One row per retained read: its id (the header up to the first tab, i.e. without
the SAM-style tags some basecallers append), the label of the ASV of its first
mapping row (or `unassigned`), the identity `100 (L − min(NM, L)) / L` with `L`
the sum of the ASV's homopolymer lengths, and the read's estimated identity.
`asvs` are the final ASVs in output order; `mappings[r]` are the rows of
`reads[r]` (ASV ids). With `upstream_quirks`, reads whose header contains a tab
are unassigned and written with the whole header, as upstream (deviation D9),
and reads sharing an id all get the assignment of the first one (upstream reads
its mapping file back by read id; deviation D13).
"""
function write_final_assignments(path::AbstractString, reads::AbstractVector{TwinRead},
                                 asvs::AbstractVector{ConsensusSequence}, mappings::AbstractVector{<:AbstractVector{ReadMapping}};
                                 upstream_quirks::Bool=false)
    final_index = Dict(c.id => i for (i, c) in enumerate(asvs))
    open(path, "w") do io
        println(io, "read_id\tasv\talignment_identity\test_read_identity")
        # The first row of a read whose ASV is still among the final ones.
        function first_final(rows)
            index = findfirst(m -> haskey(final_index, m.asv_id), rows)
            return index === nothing ? nothing : rows[index]
        end
        best = Dict{String,ReadMapping}()
        if upstream_quirks
            # As upstream reads its mapping file back: by read id, the first read
            # with an id wins; a tab in the header makes the ASV field unreadable (D9).
            for (read, rows) in zip(reads, mappings)
                (occursin('\t', read.id) || haskey(best, read.id)) && continue
                m = first_final(rows)
                m === nothing || (best[read.id] = m)
            end
        end
        for (read, rows) in zip(reads, mappings)
            estimated = read.estimated_identity === nothing ? "NA" : @sprintf("%.4f", read.estimated_identity)
            mapping = upstream_quirks ? get(best, read.id, nothing) : first_final(rows)
            name = upstream_quirks ? read.id : first(split(read.id, '\t'))
            if mapping === nothing
                println(io, name, "\tunassigned\tNA\t", estimated)
            else
                i = final_index[mapping.asv_id]; asv = asvs[i]
                len = sum(Int, asv.hp_lengths)
                identity = len > 0 ? (len - min(mapping.nm, len)) / len * 100.0 : 0.0
                println(io, name, "\tfinal_asv_", i - 1, "_depth_", _depth_field(asv), '\t', @sprintf("%.4f", identity), '\t', estimated)
            end
        end
    end
    return path
end

"""
    savont_asv(inputs, config=SavontConfig(); output_dir=nothing) -> NamedTuple

Run the whole `savont asv` pipeline on the FASTQ/FASTA files `inputs`: markers
(stage 1), retained reads, primary and polymorphic-marker clustering (stages
2–3), initial consensus, pileups and polishing (stage 4), merging (stage 5),
chimera detection (stage 6, unless `skip_chimera_detection`), EM quantification
(stage 7) and, with `pooled_samples` and several inputs, per-sample depths. The
ASVs are sorted by decreasing depth, then sequence, then cluster index, as
upstream.

With `output_dir`, writes `final_asvs.fasta`, `feature-table.tsv` and
`final_assignments.tsv` there. Returns `(asvs, reads, mappings, low_quality,
chimeras, low_polymorphism)`: the final ASVs, the retained reads, each read's
assignment rows, the consensuses set aside as low quality, the chimera calls
(indices into the stage-5 consensuses) and whether low-polymorphism mode was in
effect.

Upstream ends the run after stage 5, without quantification or output files,
when chimera detection is skipped (deviation D8); with `upstream_quirks` this
function does the same and returns `asvs = nothing`.
"""
function savont_asv(inputs::AbstractVector{<:AbstractString}, config::SavontConfig=SavontConfig();
                    output_dir::Union{Nothing,AbstractString}=nothing)
    validate(config)
    paths = collect(String, inputs)
    info, blockmers = discover_markers(paths, config)
    reads, detected_low = build_twin_reads(paths, info, blockmers, config)
    low_polymorphism = config.low_polymorphism || detected_low
    detected_low && !config.low_polymorphism &&
        @warn "Auto-enabling low-polymorphism mode: >75% of reads have no SNPmers."
    primary = cluster_reads_by_kmers(reads, config)
    clusters = cluster_reads_by_snpmers(reads, primary, config; low_polymorphism)
    consensuses = initial_consensus(reads, clusters, config)
    kept, low_quality, _ = polish_consensuses!(reads, consensuses, config)
    merged = merge_similar_consensuses(kept, config)
    chimeras = ChimeraInfo[]
    if config.skip_chimera_detection
        config.upstream_quirks && return (asvs=nothing, reads, mappings=nothing, low_quality, chimeras, low_polymorphism)
        remaining = merged
    else
        chimeras = detect_chimeras(merged, config)
        remaining = filter_chimeras(merged, chimeras)
    end
    asvs, mappings = refine_depths!(reads, remaining, info, config; low_polymorphism)
    sort!(asvs; by=c -> (-c.depth, c.sequence, c.id))
    if config.pooled_samples && length(paths) > 1
        depths = per_sample_depths(reads, asvs, info, config, length(paths); low_polymorphism)
        for (i, c) in enumerate(asvs)
            c.per_sample_depths = depths[i, :]
        end
    end
    if output_dir !== nothing
        mkpath(output_dir)
        write_asv_fasta(joinpath(output_dir, ASV_FILE), asvs, "final")
        # Rust `Path::file_stem`: the name without its last extension.
        names = [first(splitext(basename(p))) for p in paths]
        write_feature_table(joinpath(output_dir, "feature-table.tsv"), asvs,
                            config.pooled_samples && length(paths) > 1 ? names : names[1:1])
        write_final_assignments(joinpath(output_dir, "final_assignments.tsv"), reads, asvs, mappings;
                                upstream_quirks=config.upstream_quirks)
    end
    return (; asvs, reads, mappings, low_quality, chimeras, low_polymorphism)
end
