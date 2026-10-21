# Phase 3C, second part: primary k-mer clustering (upstream
# `asv_cluster::cluster_reads_by_kmers`, MinHash LSH path).

const PRIMARY_TOP_CANDIDATES = 10

"""
    cluster_reads_by_kmers(reads::Vector{TwinRead}, config::SavontConfig) -> Vector{Vector{Int}}

Greedy primary clustering (upstream stage 2). Reads are visited in order (as
returned by [`build_twin_reads`](@ref)); each is compared with the existing
representatives that share at least one LSH bucket signature. Candidates are
ranked by decreasing number of shared buckets, then decreasing index; those with
the maximal count, and at least the first ten, are verified with

    similarity = (shared distinct minimizers / max(distinct read minimizers,
                  representative minimizers)) ^ (1 / k)

The read joins the first candidate with the highest similarity when it exceeds
`primary_clustering_threshold`, otherwise it becomes a new representative.

Returns clusters as sorted vectors of 1-based read indices, ordered by decreasing
size and then by first member, without clusters smaller than `min_cluster_size`.
"""
function cluster_reads_by_kmers(reads::AbstractVector{TwinRead}, config::SavontConfig)
    k = config.kmer_size
    threshold = config.primary_clustering_threshold
    bucket_index = [Dict{UInt64,Vector{Int}}() for _ in 1:LSH_NUM_TABLES]
    representative_sets = Dict{Int,Set{UInt64}}()
    assignment = zeros(Int, length(reads))
    hits = Dict{Int,Int}()
    for (read_index, read) in enumerate(reads)
        empty!(hits)
        for (table, signature) in enumerate(read.lsh_signatures)
            signature === nothing && continue
            for candidate in get(bucket_index[table], signature, ())
                hits[candidate] = get(hits, candidate, 0) + 1
            end
        end
        best = 0
        if !isempty(hits)
            candidates = sort!(collect(hits); by=x -> (x[2], x[1]), rev=true)
            maximum_hits = candidates[1][2]
            read_set = Set(read.minimizer_kmers)
            best_similarity = 0.0
            for (checked, (candidate, candidate_hits)) in enumerate(candidates)
                candidate_hits == maximum_hits || checked <= PRIMARY_TOP_CANDIDATES || break
                representative_set = representative_sets[candidate]
                shared = count(in(representative_set), read_set)
                ratio = shared / max(length(read_set), length(reads[candidate].minimizer_kmers))
                similarity = ratio^(1.0 / k)
                if similarity > best_similarity
                    best_similarity = similarity
                    best = candidate
                end
            end
            best_similarity > threshold || (best = 0)
        end
        if best == 0
            for (table, signature) in enumerate(read.lsh_signatures)
                signature === nothing || push!(get!(Vector{Int}, bucket_index[table], signature), read_index)
            end
            representative_sets[read_index] = Set(read.minimizer_kmers)
            best = read_index
        end
        assignment[read_index] = best
    end
    members = Dict{Int,Vector{Int}}()
    for (read_index, representative) in enumerate(assignment)
        push!(get!(Vector{Int}, members, representative), read_index)
    end
    clusters = sort!(collect(values(members)); by=cluster -> (-length(cluster), first(cluster)))
    return filter!(cluster -> length(cluster) >= config.min_cluster_size, clusters)
end
