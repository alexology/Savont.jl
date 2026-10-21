# Partial order alignment, ported from spoars 0.1.3 (a bit-exact Rust port of
# spoa 4.1.5) for the configuration upstream savont uses: overlap alignment with a
# linear gap penalty, optionally banded. Node and edge identifiers are 1-based.

"""
    POAScoring(match, mismatch, gap)

Linear-gap scoring for [`poa_align`](@ref): `match` for identical bases,
`mismatch` otherwise and `gap` per inserted or deleted base (`gap <= 0`).
"""
struct POAScoring
    match::Int
    mismatch::Int
    gap::Int
    function POAScoring(match::Integer, mismatch::Integer, gap::Integer)
        gap <= 0 || throw(ArgumentError("gap penalty must be non-positive"))
        return new(match, mismatch, gap)
    end
end

"""
    POABand(base, fraction)

abPOA-style adaptive band (spoars `BandConfig`): each graph row computes only the
query columns within `base + round(fraction * length)` of its anchor. `fraction`
is stored as a `Float32`, as upstream.
"""
struct POABand
    base::Int
    fraction::Float32
end
_band_width(band::POABand, query_length::Integer) =
    band.base + max(Int(round(Float64(band.fraction) * query_length, RoundNearestTiesAway)), 0)

mutable struct POAGraph
    codes::Vector{Int}                 # 0-based symbol code of each node
    inedges::Vector{Vector{Int}}
    outedges::Vector{Vector{Int}}
    aligned::Vector{Vector{Int}}
    edge_tail::Vector{Int}
    edge_head::Vector{Int}
    edge_labels::Vector{Vector{Int}}
    edge_weight::Vector{Int64}
    coder::Vector{Int}                 # byte + 1 => code, or -1
    decoder::Vector{Int}               # code + 1 => byte
    sequences::Vector{Int}             # first node of each added sequence
    rank_to_node::Vector{Int}
end

"""
    POAGraph()

An empty partial order alignment graph (spoa `Graph`).
"""
POAGraph() = POAGraph(Int[], Vector{Int}[], Vector{Int}[], Vector{Int}[], Int[], Int[], Vector{Int}[], Int64[],
                      fill(-1, 256), Int[], Int[], Int[])

Base.length(graph::POAGraph) = length(graph.codes)

function _add_node!(graph::POAGraph, code::Int)
    push!(graph.codes, code); push!(graph.inedges, Int[]); push!(graph.outedges, Int[]); push!(graph.aligned, Int[])
    return length(graph.codes)
end

function _add_edge!(graph::POAGraph, tail::Int, head::Int, weight::Int64, label::Int)
    for edge in graph.outedges[tail]
        if graph.edge_head[edge] == head
            push!(graph.edge_labels[edge], label)
            graph.edge_weight[edge] += weight
            return nothing
        end
    end
    push!(graph.edge_tail, tail); push!(graph.edge_head, head)
    push!(graph.edge_labels, [label]); push!(graph.edge_weight, weight)
    edge = length(graph.edge_tail)
    push!(graph.outedges[tail], edge); push!(graph.inedges[head], edge)
    return nothing
end

# Positions `first:last` (1-based) of `sequence` as a new chain; returns its first node or 0.
function _add_sequence!(graph::POAGraph, sequence, weights, first::Int, last::Int)
    first > last && return 0
    label = length(graph.sequences)
    start = length(graph.codes) + 1
    previous = 0
    for i in first:last
        current = _add_node!(graph, graph.coder[sequence[i] + 1])
        previous == 0 || _add_edge!(graph, previous, current, Int64(weights[i - 1]) + Int64(weights[i]), label)
        previous = current
    end
    return start
end

"""
    poa_add!(graph, alignment, sequence, weights)

Add `sequence` to `graph` along `alignment` (spoa `Graph::AddAlignment`).
`alignment` holds `(node, position)` pairs with 0 for a gap on either side, as
returned by [`poa_align`](@ref); an empty alignment adds the sequence as a new
chain. Edge weights are the sums of the `weights` of consecutive bases.
"""
function poa_add!(graph::POAGraph, alignment::AbstractVector{Tuple{Int,Int}}, sequence::AbstractVector{UInt8},
                  weights::AbstractVector{<:Integer})
    isempty(sequence) && return graph
    length(sequence) == length(weights) || throw(ArgumentError("sequence and weights lengths differ"))
    for base in sequence
        if graph.coder[base + 1] == -1
            graph.coder[base + 1] = length(graph.decoder)
            push!(graph.decoder, Int(base))
        end
    end
    if isempty(alignment)
        start = _add_sequence!(graph, sequence, weights, 1, length(sequence))
        start == 0 || push!(graph.sequences, start)
        _topological_sort!(graph)
        return graph
    end
    positions = [position for (_, position) in alignment if position != 0]
    isempty(positions) && throw(ArgumentError("alignment contains no sequence position"))
    front = first(positions); back = last(positions)
    begin_node = _add_sequence!(graph, sequence, weights, 1, front - 1)
    previous = begin_node == 0 ? 0 : length(graph.codes)
    suffix = _add_sequence!(graph, sequence, weights, back + 1, length(sequence))
    label = length(graph.sequences)
    for (node, position) in alignment
        position == 0 && continue
        code = graph.coder[sequence[position] + 1]
        current = if node == 0
            _add_node!(graph, code)
        elseif graph.codes[node] == code
            node
        else
            matched = 0
            for other in graph.aligned[node]
                if graph.codes[other] == code
                    matched = other
                    break
                end
            end
            if matched == 0
                matched = _add_node!(graph, code)
                for other in copy(graph.aligned[node])
                    push!(graph.aligned[other], matched); push!(graph.aligned[matched], other)
                end
                push!(graph.aligned[node], matched); push!(graph.aligned[matched], node)
            end
            matched
        end
        begin_node == 0 && (begin_node = current)
        previous == 0 || _add_edge!(graph, previous, current, Int64(weights[position - 1]) + Int64(weights[position]), label)
        previous = current
    end
    if suffix != 0
        _add_edge!(graph, previous, suffix, Int64(weights[back]) + Int64(weights[back + 1]), label)
    end
    begin_node == 0 || push!(graph.sequences, begin_node)
    _topological_sort!(graph)
    return graph
end

function _topological_sort!(graph::POAGraph)
    empty!(graph.rank_to_node)
    n = length(graph.codes)
    marks = zeros(UInt8, n); ignored = falses(n)
    stack = Int[]
    for start in 1:n
        marks[start] != 0 && continue
        push!(stack, start)
        while !isempty(stack)
            current = stack[end]
            valid = true
            if marks[current] != 2
                for edge in graph.inedges[current]
                    tail = graph.edge_tail[edge]
                    if marks[tail] != 2
                        push!(stack, tail)
                        valid = false
                    end
                end
                if !ignored[current]
                    for other in graph.aligned[current]
                        if marks[other] != 2
                            push!(stack, other)
                            ignored[other] = true
                            valid = false
                        end
                    end
                end
                if valid
                    marks[current] = 2
                    if !ignored[current]
                        push!(graph.rank_to_node, current)
                        append!(graph.rank_to_node, graph.aligned[current])
                    end
                else
                    marks[current] = 1
                end
            end
            valid && pop!(stack)
        end
    end
    return graph
end

const _POA_NEG_INF = Int(typemin(Int32)) + 1024

# Scalar cell readers for the backtrack (spoars `CellRead`).
struct _ExactMatrix
    h::Matrix{Int}          # (seq_len + 1) × (nodes + 1)
end
@inline _cell(m::_ExactMatrix, i, j) = m.h[j + 1, i + 1]

struct _StripedMatrix{T}
    h::Matrix{T}            # padded columns × (nodes + 1); row 0 and column 0 are the boundary
    gap::Int
end
@inline _cell(m::_StripedMatrix, i, j) = i == 0 ? j * m.gap : j == 0 ? 0 : Int(m.h[j, i + 1])

_profile_score(graph, scoring, code, base) = graph.decoder[code + 1] == base ? scoring.match : scoring.mismatch

function _backtrack(graph::POAGraph, sequence, scoring::POAScoring, h, rank_of, max_i::Int, max_j::Int)
    alignment = Tuple{Int,Int}[]
    (max_i == 0 && max_j == 0) && return alignment
    g = scoring.gap
    pred_row(edge) = rank_of[graph.edge_tail[edge]] + 1
    i = max_i; j = max_j
    while !(i == 0 || j == 0)
        h_ij = _cell(h, i, j)
        node = graph.rank_to_node[i]
        inedges = graph.inedges[node]
        previous_i = 0; previous_j = 0; found = false
        match_score = _profile_score(graph, scoring, graph.codes[node], sequence[j])
        first_pred = isempty(inedges) ? 0 : pred_row(inedges[1])
        if h_ij == _cell(h, first_pred, j - 1) + match_score
            previous_i, previous_j, found = first_pred, j - 1, true
        else
            for p in 2:length(inedges)
                pred = pred_row(inedges[p])
                if h_ij == _cell(h, pred, j - 1) + match_score
                    previous_i, previous_j, found = pred, j - 1, true
                    break
                end
            end
        end
        if !found
            if h_ij == _cell(h, first_pred, j) + g
                previous_i, previous_j, found = first_pred, j, true
            else
                for p in 2:length(inedges)
                    pred = pred_row(inedges[p])
                    if h_ij == _cell(h, pred, j) + g
                        previous_i, previous_j, found = pred, j, true
                        break
                    end
                end
            end
        end
        if !found && h_ij == _cell(h, i, j - 1) + g
            previous_i, previous_j = i, j - 1
        end
        push!(alignment, (i == previous_i ? 0 : node, j == previous_j ? 0 : j))
        i, j = previous_i, previous_j
    end
    return reverse!(alignment)
end

function _worst_case_score(scoring::POAScoring, i::Int, j::Int)
    gap_score(len) = len == 0 ? 0 : min(scoring.gap * len, 0)
    return min(-(scoring.match * min(i, j) + gap_score(abs(i - j))), gap_score(i) + gap_score(j))
end

# Exact overlap alignment (spoars `SisdEngine::linear`).
function _poa_align_exact(graph::POAGraph, sequence, scoring::POAScoring, rank_of)
    n = length(sequence); rows = length(graph) + 1
    g = scoring.gap
    h = Matrix{Int}(undef, n + 1, rows)
    h[1, 1] = 0
    for j in 1:n
        h[j + 1, 1] = j * g
    end
    for r in 2:rows
        h[1, r] = 0
    end
    max_score = _POA_NEG_INF; max_i = 0; max_j = 0
    for node in graph.rank_to_node
        i = rank_of[node] + 1
        code = graph.codes[node]
        inedges = graph.inedges[node]
        pred = isempty(inedges) ? 0 : rank_of[graph.edge_tail[inedges[1]]] + 1
        for j in 1:n
            h[j + 1, i + 1] = max(h[j, pred + 1] + _profile_score(graph, scoring, code, sequence[j]), h[j + 1, pred + 1] + g)
        end
        for p in 2:length(inedges)
            pred = rank_of[graph.edge_tail[inedges[p]]] + 1
            for j in 1:n
                m = h[j, pred + 1] + _profile_score(graph, scoring, code, sequence[j])
                h[j + 1, i + 1] = max(m, max(h[j + 1, i + 1], h[j + 1, pred + 1] + g))
            end
        end
        sink = isempty(graph.outedges[node])
        for j in 1:n
            value = max(h[j, i + 1] + g, h[j + 1, i + 1])
            h[j + 1, i + 1] = value
            if sink && max_score < value
                max_score, max_i, max_j = value, i, j
            end
        end
    end
    return _backtrack(graph, sequence, scoring, _ExactMatrix(h), rank_of, max_i, max_j), max_score
end

# Remaining heaviest-path length per rank (spoars `band::remaining_path`).
function _remaining_path(graph::POAGraph, rank_of)
    r = zeros(Int, length(graph.rank_to_node))
    for node in Iterators.reverse(graph.rank_to_node)
        best_weight = Int64(0); best_rank = -1
        for edge in graph.outedges[node]
            successor = rank_of[graph.edge_head[edge]]
            if best_rank == -1 || graph.edge_weight[edge] > best_weight ||
               (graph.edge_weight[edge] == best_weight && successor < best_rank)
                best_weight, best_rank = graph.edge_weight[edge], successor
            end
        end
        r[rank_of[node] + 1] = best_rank == -1 ? 0 : 1 + r[best_rank + 1]
    end
    return r
end

@inline _sat(::Type{T}, x::Int) where {T} = T(clamp(x, Int(typemin(T)), Int(typemax(T))))
@inline _adds(a::T, b::T) where {T} = _sat(T, Int(a) + Int(b))

# Banded overlap alignment emulating spoars' SIMD kernels lane by lane: striped
# segments of `lanes` query columns, element type `T` (Int16 or Int32) with the
# kernels' saturating and wrapping arithmetic, and spoars' abPOA band geometry.
function _poa_align_banded(::Type{T}, lanes::Int, graph::POAGraph, sequence, scoring::POAScoring,
                           rank_of, band::POABand) where {T}
    n = length(sequence)
    segments = cld(n, lanes); padded = segments * lanes
    rows = length(graph) + 1
    neg = T(Int(typemin(T)) + 1024)
    g = T(scoring.gap)
    padding = -max(abs(scoring.match), abs(scoring.mismatch), abs(scoring.gap))
    profile = Matrix{T}(undef, padded, length(graph.decoder))
    for code in 0:(length(graph.decoder) - 1), position in 1:padded
        profile[position, code + 1] = position <= n ? T(_profile_score(graph, scoring, code, sequence[position])) : T(padding)
    end
    h = fill(neg, padded, rows)
    for position in 1:n
        h[position, 1] = T(position * scoring.gap)
    end
    remaining = _remaining_path(graph, rank_of)
    best_col = zeros(Int, length(graph.rank_to_node))
    width = _band_width(band, n)
    steps = trailing_zeros(lanes)
    max_score = _POA_NEG_INF; max_i = 0
    floor_value = Int(neg)
    hv = Vector{T}(undef, lanes)
    for node in graph.rank_to_node
        rank = rank_of[node]; i = rank + 1
        inedges = graph.inedges[node]
        code = graph.codes[node]
        # Band window in segments (spoars `row_band`, `node_window`, `segment_range`).
        anchor = max(n - remaining[rank + 1], 0)
        if isempty(inedges)
            low = high = anchor
        else
            low = typemax(Int); high = 0
            for edge in inedges
                column = best_col[rank_of[graph.edge_tail[edge]] + 1]
                low = min(low, column); high = max(high, column)
        end
        end
        window_begin = max(min(low, anchor) - width, 0)
        window_end = min(max(high, anchor) + width + 1, n)
        first_segment = min(window_begin ÷ lanes, segments - 1)
        last_segment = min(cld(window_end, lanes), segments)
        last_segment <= first_segment && (last_segment = min(first_segment + 1, segments))
        first_column = first_segment * lanes + 1
        last_column = last_segment * lanes
        # Diagonal and vertical moves, first predecessor then the others.
        for (p, edge) in enumerate(isempty(inedges) ? [0] : inedges)
            pred = edge == 0 ? 0 : rank_of[graph.edge_tail[edge]] + 1
            carry = first_segment == 0 ? T(0) : h[first_column - 1, pred + 1]
            for position in first_column:last_column
                diagonal = position == first_column ? carry : h[position - 1, pred + 1]
                candidate = max(_adds(diagonal, profile[position, code + 1]), _adds(h[position, pred + 1], g))
                h[position, i + 1] = p == 1 ? candidate : max(h[position, i + 1], candidate)
            end
        end
        # Horizontal gaps: per segment, lane-0 carry then the prefix-max ladder.
        carry = first_segment == 0 ? g : neg      # first_column(i) + g, or closed edge
        score = fill(neg, lanes)
        for segment in first_segment:(last_segment - 1)
            offset = segment * lanes
            for k in 1:lanes
                hv[k] = max(h[offset + k, i + 1], k == 1 ? carry : neg)
            end
            for step in 0:(steps - 1)
                shift = 1 << step
                penalty = T(Int(g) << step)
                for k in lanes:-1:1
                    shifted = k > shift ? hv[k - shift] + penalty : neg     # wrapping add
                    hv[k] = max(hv[k], shifted)
                end
            end
            carry = _adds(hv[lanes], g)
            for k in 1:lanes
                h[offset + k, i + 1] = hv[k]
                score[k] = max(score[k], hv[k])
            end
        end
        row_best = Int(maximum(score))
        column = findfirst(==(row_best), @view h[first_column:last_column, i + 1])
        best_col[rank + 1] = max((first_column - 1) + (column === nothing ? -1 : column - 1), 0)
        if isempty(graph.outedges[node])
            row_score = row_best <= floor_value ? _POA_NEG_INF : row_best
            if max_score < row_score
                max_score, max_i = row_score, i
            end
        end
    end
    max_i == 0 && return Tuple{Int,Int}[], max_score
    max_score == _POA_NEG_INF && return Tuple{Int,Int}[], _POA_NEG_INF
    column = findfirst(x -> Int(x) == max_score, @view h[:, max_i + 1])
    column === nothing && return Tuple{Int,Int}[], max_score
    return _backtrack(graph, sequence, scoring, _StripedMatrix{T}(h, scoring.gap), rank_of, max_i, column), max_score
end

"""
    poa_align(graph, sequence, scoring; band=nothing, isa=:avx2) -> (alignment, score)

Overlap alignment of `sequence` to `graph` with a linear gap penalty (spoars
`SimdEngine`). Without `band` the alignment is exact (spoa's). With a band, the
result reproduces spoars' banded SIMD kernels, whose computed region depends on
the vector width: `isa = :avx2` (16 lanes for 16-bit scores, 8 for 32-bit) or
`:sse41` (8 and 4; also NEON). Scores switch from 16 to 32 bits when the
sequence and graph are long enough to risk overflow; when even 32 bits could
overflow, the exact alignment is used. An empty alignment means the sequence
should be added as a new chain.
"""
function poa_align(graph::POAGraph, sequence::AbstractVector{UInt8}, scoring::POAScoring;
                   band::Union{Nothing,POABand}=nothing, isa::Symbol=:avx2)
    (isempty(sequence) || length(graph) == 0) && return Tuple{Int,Int}[], 0
    isa in (:avx2, :sse41) || throw(ArgumentError("isa must be :avx2 or :sse41"))
    rank_of = zeros(Int, length(graph))
    for (rank, node) in enumerate(graph.rank_to_node)
        rank_of[node] = rank - 1
    end
    worst = _worst_case_score(scoring, length(sequence) + 8, length(graph))
    if band === nothing || worst < Int(typemin(Int32)) + 1024
        _worst_case_score(scoring, length(sequence), length(graph)) >= _POA_NEG_INF ||
            throw(ArgumentError("possible overflow"))
        return _poa_align_exact(graph, sequence, scoring, rank_of)
    end
    wide = isa === :avx2
    if worst < Int(typemin(Int16)) + 1024
        return _poa_align_banded(Int32, wide ? 8 : 4, graph, sequence, scoring, rank_of, band)
    end
    return _poa_align_banded(Int16, wide ? 16 : 8, graph, sequence, scoring, rank_of, band)
end

"""
    poa_consensus(graph) -> Vector{UInt8}

Heaviest-bundle consensus of `graph` (spoa `Graph::GenerateConsensus`).
"""
function poa_consensus(graph::POAGraph)
    isempty(graph.rank_to_node) && return UInt8[]
    n = length(graph)
    predecessors = zeros(Int, n); scores = fill(Int64(-1), n)
    function relax!(node, skip_unscored)
        for edge in graph.inedges[node]
            tail = graph.edge_tail[edge]
            skip_unscored && scores[tail] == -1 && continue
            weight = graph.edge_weight[edge]
            if scores[node] < weight || (scores[node] == weight && scores[predecessors[node]] <= scores[tail])
                scores[node] = weight
                predecessors[node] = tail
            end
        end
        predecessors[node] == 0 || (scores[node] += scores[predecessors[node]])
    end
    best = 0
    for node in graph.rank_to_node
        relax!(node, false)
        (best == 0 || scores[best] < scores[node]) && (best = node)
    end
    if !isempty(graph.outedges[best])
        rank_of = zeros(Int, n)
        for (rank, node) in enumerate(graph.rank_to_node)
            rank_of[node] = rank
        end
        while !isempty(graph.outedges[best])
            start = best
            for edge in graph.outedges[start], in_edge in graph.inedges[graph.edge_head[edge]]
                tail = graph.edge_tail[in_edge]
                tail == start || (scores[tail] = -1)
            end
            best = 0
            for rank in (rank_of[start] + 1):length(graph.rank_to_node)
                node = graph.rank_to_node[rank]
                scores[node] = -1; predecessors[node] = 0
                relax!(node, true)
                (best == 0 || scores[best] < scores[node]) && (best = node)
            end
        end
    end
    path = Int[]
    node = best
    while predecessors[node] != 0
        push!(path, node)
        node = predecessors[node]
    end
    push!(path, node)
    return UInt8[graph.decoder[graph.codes[node] + 1] for node in Iterators.reverse(path)]
end
