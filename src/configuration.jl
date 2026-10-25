Base.@kwdef struct SavontConfig
    min_read_length::Int = 1100
    max_read_length::Int = 2000
    quality_value_cutoff::Float64 = 98.0
    minimum_base_quality::UInt8 = 25
    single_strand::Bool = false
    min_cluster_size::Int = 12
    bloom_filter_size::Float64 = 0.0
    n_depth_cutoff::Int = 250
    use_hpc::Bool = false
    mask_low_quality::Bool = false
    posterior_threshold_ln::Float64 = 30.0
    max_iterations_recluster::Int = 10
    aggressive_bloom::Bool = false
    skip_chimera_detection::Bool = false
    no_snpmers::Bool = false
    low_polymorphism::Bool = false
    kmer_size::Int = 17
    blockmer_length::Int = 3
    use_blockmers::Bool = false
    chimera_allowable_errors::Int = 1
    chimera_detect_length::Union{Nothing,Int} = nothing
    phase_heterogeneous::Bool = false
    no_band::Bool = false
    primary_clustering_threshold::Float64 = 0.95
    sampling_rate::Int = 11
    pooled_samples::Bool = false
    # Vector width emulated by the banded partial order alignment (:avx2 or :sse41).
    poa_isa::Symbol = :avx2
    # Reproduce upstream bugs listed in docs/deviations.md (for oracle comparisons).
    upstream_quirks::Bool = false
end

function validate(config::SavontConfig)
    config.min_read_length >= 0 || throw(ArgumentError("min_read_length must be nonnegative"))
    config.max_read_length >= config.min_read_length || throw(ArgumentError("max_read_length must be at least min_read_length"))
    0.0 <= config.quality_value_cutoff <= 100.0 || throw(ArgumentError("quality_value_cutoff must be in [0, 100]"))
    config.min_cluster_size > 0 || throw(ArgumentError("min_cluster_size must be positive"))
    isodd(config.kmer_size) && 0 < config.kmer_size < 24 || throw(ArgumentError("kmer_size must be odd and smaller than 24"))
    config.kmer_size + config.blockmer_length <= 31 || throw(ArgumentError("kmer_size + blockmer_length must be at most 31"))
    0 < config.sampling_rate <= config.kmer_size || throw(ArgumentError("sampling_rate must be in 1:kmer_size"))
    # An even sampling rate has no central s-mer, so the two strands of a read get
    # different syncmers (deviation D11); upstream accepts it.
    isodd(config.sampling_rate) || config.upstream_quirks ||
        throw(ArgumentError("sampling_rate must be odd (an even value separates the two strands)"))
    0.0 <= config.primary_clustering_threshold <= 1.0 || throw(ArgumentError("primary_clustering_threshold must be in [0, 1]"))
    return config
end

function apply_preset(config::SavontConfig, preset::Symbol)
    values = (; (field => getfield(config, field) for field in fieldnames(SavontConfig))...)
    updated = if preset === :fl_16s
        values
    elseif preset === :hifi
        merge(values, (; min_cluster_size=4))
    elseif preset === :rrna_operon
        merge(values, (; min_read_length=3500, max_read_length=5000))
    else
        throw(ArgumentError("unknown preset: $preset"))
    end
    return validate(SavontConfig(; updated...))
end

# Threaded map keeping the order of `items` (one task per item). Every task is
# waited for before the first failure (in item order) is rethrown, so no task is
# still running when the caller returns or throws.
function _tmap(f, items)
    tasks = [Threads.@spawn f(item) for item in items]
    failure = nothing
    for task in tasks
        try
            wait(task)
        catch e
            failure === nothing && (failure = e)
        end
    end
    failure === nothing || throw(failure)
    return map(fetch, tasks)
end
