# Store the upstream oracle outputs in test/reference.
#
#   SAVONT_RUST_ORACLE=oracle/rust_upstream/target/release/savont-upstream-oracle.exe \
#       julia --project=. tools/make_reference.jl
#
# Build the oracle first (see oracle/README.md).

using CodecZlib

include(joinpath(@__DIR__, "..", "test", "oracle_cases.jl"))

function store(command, name)
    text = read(command, String)
    isempty(text) && error("empty oracle output for $name")
    write(joinpath(REFERENCE_DIRECTORY, name), transcode(GzipCompressor, Vector{UInt8}(text)))
    println(name, ": ", count(==('\n'), text), " lines")
end

function main()
    oracle = get(ENV, "SAVONT_RUST_ORACLE", "")
    isempty(oracle) && error("set SAVONT_RUST_ORACLE to the compiled oracle executable")
    mkpath(REFERENCE_DIRECTORY)
    for (name, input, db, arguments, _) in CLASSIFY_RUNS
        store(classify_command(oracle, input, db, arguments), classify_reference_name(name))
    end
    for (name, input, db, arguments, _) in SINTAX_RUNS
        store(sintax_command(oracle, input, db, arguments), sintax_reference_name(name))
    end
    for (name, inputs, arguments, _) in EXPORT_RUNS
        store(export_command(oracle, inputs, arguments), export_reference_name(name))
    end
    # `julia tools/make_reference.jl taxonomy` stores only the taxonomy cases.
    "taxonomy" in ARGS && return
    for input in PRIMITIVE_INPUTS, p in PRIMITIVE_PARAMETERS
        store(primitive_command(oracle, input, p), primitive_reference_name(input, p))
    end
    store(stats_command(oracle), STATS_REFERENCE)
    for (name, arguments, _) in MARKER_RUNS
        store(marker_command(oracle, arguments), marker_reference_name(name))
    end
    for (input, mode, isa) in POA_RUNS
        store(poa_command(oracle, input, mode, isa), poa_reference_name(input, mode, isa))
    end
    for (name, inputs, arguments, _) in ASV_RUNS
        store(asv_command(oracle, inputs, arguments), asv_reference_name(name))
    end
    for (name, inputs, arguments, _) in CLUSTER_RUNS
        store(cluster_command(oracle, inputs, arguments), cluster_reference_name(name))
        store(consensus_command(oracle, inputs, arguments), consensus_reference_name(name))
        store(polish_command(oracle, inputs, arguments), polish_reference_name(name))
        store(merge_command(oracle, inputs, arguments), merge_reference_name(name))
        store(chimera_command(oracle, inputs, arguments), chimera_reference_name(name))
    end
end

main()
