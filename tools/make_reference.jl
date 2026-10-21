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
    for input in PRIMITIVE_INPUTS, p in PRIMITIVE_PARAMETERS
        store(primitive_command(oracle, input, p), primitive_reference_name(input, p))
    end
    store(stats_command(oracle), STATS_REFERENCE)
    for (name, arguments, _) in MARKER_RUNS
        store(marker_command(oracle, arguments), marker_reference_name(name))
    end
    for (name, inputs, arguments, _) in CLUSTER_RUNS
        store(cluster_command(oracle, inputs, arguments), cluster_reference_name(name))
    end
end

main()
