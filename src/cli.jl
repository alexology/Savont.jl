# Command-line interface with the subcommands and options of savont 0.8.0
# (upstream `cli.rs` and `main.rs`).

const UPSTREAM_VERSION = "0.8.0"

# One option: long name => (short letter or nothing, kind, default, help).
# Kinds: :flag, :int, :float, :string, :strings (one or more values).
const _Option = Tuple{Union{Nothing,Char},Symbol,Any,String}

const ASV_OPTIONS = Dict{String,_Option}(
    "output-dir" => ('o', :string, "savont-out", "Output directory for results (created if it does not exist)"),
    "threads" => ('t', :int, nothing, "Accepted for compatibility; the Julia threads are set with `julia -t N`"),
    "fl-16s" => (nothing, :flag, false, "16S rRNA full length (~1500 bp) amplicon preset (default; does nothing)"),
    "hifi" => (nothing, :flag, false, "PacBio HiFi amplicon preset (--min-cluster-size 4)"),
    "rrna-operon" => (nothing, :flag, false, "rRNA operon (~4000 bp) amplicon preset (--min-read-length 3500 --max-read-length 5000)"),
    "pooled-samples" => (nothing, :flag, false, "Treat each input file as a distinct sample: reads are pooled for ASV calling, depths are quantified per file"),
    "c" => ('c', :int, 11, "K-mer sampling rate: select 1 out of every C k-mers (odd)"),
    "min-read-length" => ('m', :int, 1100, "Minimum read length"),
    "max-read-length" => ('M', :int, 2000, "Maximum read length"),
    "quality-value-cutoff" => (nothing, :float, 98.0, "Minimum estimated read accuracy (%) to include in clustering"),
    "minimum-base-quality" => (nothing, :int, 25, "Minimum base quality to be considered high-quality for SNPmer detection"),
    "single-strand" => ('s', :flag, false, "Use only forward strand k-mers (for strand-specific protocols)"),
    "min-cluster-size" => (nothing, :int, 12, "Minimum number of reads required to keep a cluster (ASV)"),
    "bloom-filter-size" => ('b', :float, 0.0, "Bloom filter size in GB (0 = off; values above 0 are not implemented)"),
    "n-depth-cutoff" => ('n', :int, 250, "Minimum depth required for sequences with ambiguous bases to be included in output"),
    "use-hpc" => ('u', :flag, false, "Use homopolymer compression for clustering and consensus generation"),
    "mask-low-quality" => (nothing, :flag, false, "Mask low-quality bases in consensus sequences with N"),
    "posterior-threshold-ln" => ('p', :float, 30.0, "Negative alternate posterior probability threshold (natural log) for base consensus"),
    "max-iterations-recluster" => (nothing, :int, 10, "Maximum number of reclustering iterations"),
    "aggressive-bloom" => (nothing, :flag, false, "Accepted for compatibility (Bloom prefiltering is not implemented)"),
    "skip-chimera-detection" => (nothing, :flag, false, "Skip chimera detection (not recommended)"),
    "no-snpmers" => (nothing, :flag, false, "Disable SNPmer clustering (not recommended)"),
    "low-polymorphism" => (nothing, :flag, false, "Skip SNPmer clustering and map reads to ASVs with minimap2; enabled automatically when >75% of reads have no SNPmer"),
    "kmer-size" => ('k', :int, 17, "K-mer size for clustering (odd, < 24)"),
    "blockmer-length" => (nothing, :int, 3, "Blockmer suffix length (experimental)"),
    "use-blockmers" => (nothing, :flag, false, "Use blockmers as polymorphic markers (experimental)"),
    "chimera-allowable-errors" => (nothing, :int, 1, "Allowable errors for bi-chimeric detection"),
    "chimera-detect-length" => (nothing, :int, nothing, "Length of near-perfect segment matches for chimera detection [default: 1/10 of --min-read-length]"),
    "phase-heterogeneous" => (nothing, :flag, false, "Accepted for compatibility (unused upstream)"),
    "no-band" => (nothing, :flag, false, "Use exact (unbanded) partial order alignment for consensus generation"),
    "primary-clustering-threshold" => (nothing, :float, 0.95, "Primary k-mer clustering threshold"),
    "upstream-quirks" => (nothing, :flag, false, "Reproduce the upstream bugs listed in docs/deviations.md"),
)
const CLASSIFY_OPTIONS = Dict{String,_Option}(
    "input-dir" => ('i', :string, nothing, "Directory containing clustering results"),
    "output-dir" => ('o', :string, nothing, "Output directory [default: the input directory]"),
    "db" => ('d', :string, nothing, "Path to a savont database directory"),
    "threads" => ('t', :int, nothing, "Accepted for compatibility; the Julia threads are set with `julia -t N`"),
    "species-threshold" => (nothing, :float, 99.0, "Minimum identity (%) for species-level classification"),
    "genus-threshold" => (nothing, :float, 94.5, "Minimum identity (%) for genus-level classification"),
    "detailed-unclassified" => (nothing, :flag, false, "Write UNCLASSIFIED-(asv_header) for unassigned ranks"),
    "upstream-quirks" => (nothing, :flag, false, "Reproduce the upstream bugs listed in docs/deviations.md"),
)
const SINTAX_OPTIONS = Dict{String,_Option}(
    "input-dir" => ('i', :string, nothing, "Directory containing clustering results (with final_asvs.fasta)"),
    "output-dir" => ('o', :string, nothing, "Output directory [default: the input directory]"),
    "db" => ('d', :string, nothing, "Path to a savont database directory"),
    "threads" => ('t', :int, nothing, "Accepted for compatibility; the Julia threads are set with `julia -t N`"),
    "min-bootstrap" => (nothing, :float, 0.8, "Minimum bootstrap confidence for genus and higher ranks"),
    "n-iter" => (nothing, :int, 100, "Number of bootstrap iterations"),
    "detailed-unclassified" => (nothing, :flag, false, "Write UNCLASSIFIED-(asv_header) for unassigned ranks"),
)
const EXPORT_OPTIONS = Dict{String,_Option}(
    "input-dirs" => ('i', :strings, nothing, "Directories with savont asv (and optionally classify/sintax) output"),
    "output-dir" => ('o', :string, nothing, "Output directory for the exported results"),
    "no-fuzzy" => (nothing, :flag, false, "Keep ASVs that are identical to a longer one over up to 10 bp separate"),
    "relabel" => (nothing, :strings, nothing, "Sample names replacing the detected ones (one per sample column)"),
    "upstream-quirks" => (nothing, :flag, false, "Reproduce the upstream bugs listed in docs/deviations.md"),
)
const DOWNLOAD_OPTIONS = Dict{String,_Option}(
    "location" => (nothing, :string, nothing, "Directory under which each database is saved as a subdirectory"),
    "dbs" => (nothing, :strings, nothing, "Databases to download: $(join(DATABASE_KEYWORDS, ", "))"),
)
const SUBCOMMANDS = [
    ("asv", ASV_OPTIONS, "Turn >~98% accuracy long reads into ASVs (amplicon sequence variants)"),
    ("classify", CLASSIFY_OPTIONS, "Classify ASVs by alignment to a reference database (species/genus abundance tables)"),
    ("download", DOWNLOAD_OPTIONS, "Download reference databases (EMU, SILVA, GreenGenes2 or UNITE)"),
    ("sintax", SINTAX_OPTIONS, "Genus-level classification by k-mer bootstrap (SINTAX)"),
    ("export", EXPORT_OPTIONS, "Export ASVs to QIIME2-compatible tables, combining several runs"),
]
const REQUIRED = Dict("classify" => ["input-dir", "db"], "sintax" => ["input-dir", "db"],
                      "export" => ["input-dirs", "output-dir"], "download" => ["location", "dbs"])

struct CLIError <: Exception
    msg::String
end

function _help(command=nothing)
    io = IOBuffer()
    if command === nothing
        println(io, "Savont.jl: native Julia port of savont $UPSTREAM_VERSION, ASVs from long-read amplicons.\n")
        println(io, "Usage: julia -t <N> --project=<Savont.jl> bin/savont.jl <COMMAND> [OPTIONS]\n\nCommands:")
        foreach(((name, _, text),) -> println(io, "  ", rpad(name, 10), text), SUBCOMMANDS)
        println(io, "\nOptions:\n  -h, --help     Print help (also after a command)\n  -V, --version  Print version")
    else
        name, options, text = only(c for c in SUBCOMMANDS if c[1] == command)
        positional = name == "asv" ? " <FASTQ/FASTA>..." : ""
        println(io, text, "\n\nUsage: julia -t <N> --project=<Savont.jl> bin/savont.jl $name [OPTIONS]$positional\n\nOptions:")
        for long in sort!(collect(keys(options)))
            short, kind, default, description = options[long]
            flag = (short === nothing ? "    " : "-$short, ") * "--$long" * (kind === :flag ? "" : kind === :strings ? " <VALUE>..." : " <VALUE>")
            default === nothing || kind === :flag || (description *= " [default: $default]")
            long in get(REQUIRED, name, String[]) && (description *= " (required)")
            println(io, "  ", rpad(flag, 36), description)
        end
        println(io, "  ", rpad("-h, --help", 36), "Print help")
    end
    return String(take!(io))
end

function _parse_value(name, kind, text)
    if kind === :int
        value = occursin(r"^\+?[0-9]+$", text) ? tryparse(Int, text) : nothing
    elseif kind === :float
        value = any(isspace, text) ? nothing : tryparse(Float64, text)
    else
        value = text
    end
    value === nothing && throw(CLIError("invalid value '$text' for '--$name'"))
    return value
end

"""
    parse_command_line(args) -> (command, options, positional)

Subcommand, option values (long name => value, defaults included) and positional
arguments of a `savont` command line. Throws `Savont.CLIError` on usage errors;
`command` is `:help` or `:version` for those requests.
"""
function parse_command_line(args::AbstractVector{<:AbstractString})
    isempty(args) && throw(CLIError("a subcommand is required"))
    first_arg = String(args[1])
    first_arg in ("-h", "--help") && return (:help, Dict{String,Any}(), String[])
    first_arg in ("-V", "--version") && return (:version, Dict{String,Any}(), String[])
    index = findfirst(c -> c[1] == first_arg, SUBCOMMANDS)
    index === nothing && throw(CLIError("unrecognized subcommand '$first_arg'"))
    command, options, _ = SUBCOMMANDS[index]
    short_names = Dict(o[1] => long for (long, o) in options if o[1] !== nothing)
    values = Dict{String,Any}(); positional = String[]
    i = 2
    while i <= length(args)
        a = String(args[i])
        a in ("-h", "--help") && return (:help, Dict{String,Any}("command" => command), String[])
        if startswith(a, "--") && length(a) > 2
            name, inline = occursin('=', a) ? String.(split(a[3:end], '='; limit=2)) : (a[3:end], nothing)
        elseif startswith(a, '-') && length(a) == 2 && haskey(short_names, a[2])
            name, inline = short_names[a[2]], nothing
        elseif startswith(a, '-') && length(a) > 2 && haskey(short_names, a[2])
            name, inline = short_names[a[2]], chopprefix(a[3:end], "=")
        elseif startswith(a, '-') && a != "-"
            throw(CLIError("unexpected argument '$a' found"))
        else
            command == "asv" || throw(CLIError("unexpected argument '$a' found"))
            push!(positional, a); i += 1
            continue
        end
        haskey(options, name) || throw(CLIError("unexpected argument '$a' found"))
        haskey(values, name) && throw(CLIError("the argument '--$name' cannot be used multiple times"))
        kind = options[name][2]
        if kind === :flag
            inline === nothing || throw(CLIError("unexpected value '$inline' for '--$name' found; no more were expected"))
            values[name] = true
        elseif kind === :strings
            items = inline === nothing ? String[] : [inline]
            while i < length(args) && !(startswith(args[i + 1], '-') && args[i + 1] != "-")
                i += 1; push!(items, String(args[i]))
            end
            isempty(items) && throw(CLIError("a value is required for '--$name' but none was supplied"))
            values[name] = items
        else
            if inline === nothing
                i += 1
                i <= length(args) || throw(CLIError("a value is required for '--$name' but none was supplied"))
                inline = String(args[i])
            end
            values[name] = _parse_value(name, kind, inline)
        end
        i += 1
    end
    for name in get(REQUIRED, command, String[])
        haskey(values, name) || throw(CLIError("the following required argument was not provided: --$name"))
    end
    command == "asv" && isempty(positional) && throw(CLIError("the following required argument was not provided: <FASTQ/FASTA>..."))
    for (long, (_, kind, default, _)) in options
        haskey(values, long) || (values[long] = kind === :flag ? false : default)
    end
    return command, values, positional
end

# SavontConfig of `savont asv` options; presets are applied after the options,
# as upstream (`--rrna-operon` overrides -m/-M).
function _asv_config(values)
    config = SavontConfig(; min_read_length=values["min-read-length"], max_read_length=values["max-read-length"],
                          quality_value_cutoff=values["quality-value-cutoff"],
                          minimum_base_quality=UInt8(clamp(values["minimum-base-quality"], 0, 255)),
                          single_strand=values["single-strand"], min_cluster_size=values["min-cluster-size"],
                          bloom_filter_size=values["bloom-filter-size"], n_depth_cutoff=values["n-depth-cutoff"],
                          use_hpc=values["use-hpc"], mask_low_quality=values["mask-low-quality"],
                          posterior_threshold_ln=values["posterior-threshold-ln"],
                          max_iterations_recluster=values["max-iterations-recluster"],
                          aggressive_bloom=values["aggressive-bloom"], skip_chimera_detection=values["skip-chimera-detection"],
                          no_snpmers=values["no-snpmers"], low_polymorphism=values["low-polymorphism"],
                          kmer_size=values["kmer-size"], blockmer_length=values["blockmer-length"],
                          use_blockmers=values["use-blockmers"], chimera_allowable_errors=values["chimera-allowable-errors"],
                          chimera_detect_length=values["chimera-detect-length"], phase_heterogeneous=values["phase-heterogeneous"],
                          no_band=values["no-band"], primary_clustering_threshold=values["primary-clustering-threshold"],
                          sampling_rate=values["c"], pooled_samples=values["pooled-samples"],
                          upstream_quirks=values["upstream-quirks"])
    values["rrna-operon"] && (config = apply_preset(config, :rrna_operon))
    values["hifi"] && (config = apply_preset(config, :hifi))
    return validate(config)
end

function _run(command, values, positional, err)
    threads = get(values, "threads", nothing)
    threads === nothing || threads == Threads.nthreads() ||
        println(err, "Note: --threads $threads ignored; Savont.jl uses the $(Threads.nthreads()) Julia thread(s) of this session (julia -t N).")
    if command == "asv"
        for path in positional
            isfile(path) || throw(CLIError("Input file $path does not exist."))
        end
        output = values["output-dir"]
        ispath(output) && !isdir(output) && throw(CLIError("Output directory specified by `-o` exists and is not a directory."))
        result = savont_asv(positional, _asv_config(values); output_dir=output)
        result.asvs === nothing ? println(err, "Stopped after merging (--skip-chimera-detection with --upstream-quirks).") :
                                  println(err, "Wrote $(length(result.asvs)) ASVs to $output")
    elseif command == "classify"
        output = something(values["output-dir"], values["input-dir"])
        savont_classify(values["input-dir"], values["db"]; output_dir=output, species_threshold=values["species-threshold"],
                        genus_threshold=values["genus-threshold"], detailed_unclassified=values["detailed-unclassified"],
                        upstream_quirks=values["upstream-quirks"])
        println(err, "Wrote species_abundance.tsv, genus_abundance.tsv and asv_mappings.tsv to $output")
    elseif command == "sintax"
        output = something(values["output-dir"], values["input-dir"])
        hits = savont_sintax(values["input-dir"], values["db"]; output_dir=output, min_bootstrap=values["min-bootstrap"],
                             n_iter=values["n-iter"], detailed_unclassified=values["detailed-unclassified"])
        println(err, "SINTAX: $(length(hits)) ASVs classified; tables in $output")
    elseif command == "export"
        savont_export(values["input-dirs"], values["output-dir"]; fuzzy=!values["no-fuzzy"], relabel=values["relabel"],
                      upstream_quirks=values["upstream-quirks"])
        println(err, "Exported to $(values["output-dir"])")
    elseif command == "download"
        for keyword in values["dbs"]
            keyword in DATABASE_KEYWORDS ||
                throw(CLIError("invalid value '$keyword' for '--dbs': possible values: $(join(DATABASE_KEYWORDS, ", "))"))
        end
        for directory in download_database(values["location"], values["dbs"]...)
            println(err, "Use with: bin/savont.jl classify -d $directory")
        end
    end
end

"""
    Savont.main(args = ARGS; output = stdout, err = stderr) -> Int

Command-line entry point with the subcommands and options of savont 0.8.0
(`asv`, `classify`, `sintax`, `export`, `download`; `--help` lists them).
Returns the exit status: 0, 1 on errors, 2 on usage errors. `bin/savont.jl`
calls it.
"""
function main(args::AbstractVector{<:AbstractString}=ARGS; output::IO=stdout, err::IO=stderr)
    local command, values, positional
    try
        command, values, positional = parse_command_line(args)
    catch e
        e isa CLIError || rethrow()
        println(err, "Error: ", e.msg, "\n\nFor more information, try '--help'.")
        return 2
    end
    if command === :help
        print(output, _help(get(values, "command", nothing)))
        return 0
    elseif command === :version
        println(output, "Savont.jl ", pkgversion(@__MODULE__), " (savont ", UPSTREAM_VERSION, ")")
        return 0
    end
    try
        _run(command, values, positional, err)
    catch e
        e isa InterruptException && rethrow()
        println(err, "Error: ", e isa CLIError ? e.msg : sprint(showerror, e))
        return e isa CLIError ? 2 : 1
    end
    return 0
end
