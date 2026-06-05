# Plots.jl recipes are intentionally kept lightweight and dependency-free via RecipesBase.

function _segment_coordinates(segments)
    T = promote_type(Float64, eltype(first(segments).p))
    n = length(segments)
    x = Matrix{T}(undef, 2, n)
    y = Matrix{T}(undef, 2, n)
    z = Matrix{Int}(undef, 2, n)

    for (i, segment) in enumerate(segments)
        x[1, i] = segment.p[1]
        x[2, i] = segment.q[1]
        y[1, i] = segment.p[2]
        y[2, i] = segment.q[2]
        z[1, i] = segment.element
        z[2, i] = segment.element
    end

    return x, y, z
end

function _track_segment_coordinates(tracks)
    n = sum(track -> length(track.segments), tracks)
    iszero(n) && return Float64[], Float64[], Matrix{Int}(undef, 0, 0)

    first_segment = nothing
    for track in tracks
        if !isempty(track.segments)
            first_segment = first(track.segments)
            break
        end
    end

    T = promote_type(Float64, eltype(first_segment.p))
    x = Matrix{T}(undef, 2, n)
    y = Matrix{T}(undef, 2, n)
    z = Matrix{Int}(undef, 2, n)

    i = 0
    for track in tracks
        for segment in track.segments
            i += 1
            x[1, i] = segment.p[1]
            x[2, i] = segment.q[1]
            y[1, i] = segment.p[2]
            y[2, i] = segment.q[2]
            z[1, i] = segment.element
            z[2, i] = segment.element
        end
    end

    return x, y, z
end

@recipe function plot(track::Track)
    x = [track.p[1], track.q[1]]
    y = [track.p[2], track.q[2]]

    seriestype := :path
    linewidth --> 1.0
    legend --> false
    aspect_ratio --> :equal

    return x, y
end

@recipe function plot(tg::TrackGenerator)
    @unpack tracks_by_uid, n_total_tracks = tg

    T = promote_type(Float64, eltype(first(tracks_by_uid).p))
    x = Matrix{T}(undef, 2, n_total_tracks)
    y = Matrix{T}(undef, 2, n_total_tracks)

    for j in 1:n_total_tracks
        track = tracks_by_uid[j]
        x[1, j] = track.p[1]
        x[2, j] = track.q[1]
        y[1, j] = track.p[2]
        y[2, j] = track.q[2]
    end

    seriestype := :path
    linewidth --> 0.25
    legend --> false
    aspect_ratio --> :equal
    framestyle --> :box

    return x, y
end

@recipe function plot(segments::AbstractVector{<:Segment})
    isempty(segments) && return Float64[], Float64[]
    x, y, z = _segment_coordinates(segments)

    seriestype := :path
    linewidth --> 0.35
    legend --> false
    aspect_ratio --> :equal
    framestyle --> :box
    line_z --> z

    return x, y
end

@recipe function plot(tracks::AbstractVector{<:Track})
    x, y, z = _track_segment_coordinates(tracks)

    seriestype := :path
    linewidth --> 0.35
    legend --> false
    aspect_ratio --> :equal
    framestyle --> :box
    line_z --> z

    return x, y
end

@recipe function plot(mesh::Mesh)
    @unpack cell_nodes, model = mesh

    grid = get_grid(model)
    nodes = get_node_coordinates(grid)

    n_cells = length(cell_nodes)
    max_nodes = maximum(length, cell_nodes)

    x = Matrix{Float64}(undef, max_nodes + 1, n_cells)
    y = Matrix{Float64}(undef, max_nodes + 1, n_cells)
    fill!(x, NaN)
    fill!(y, NaN)

    for (i, node_ids) in enumerate(cell_nodes)
        ids = ordered_node_ids(mesh, node_ids)
        for (j, node_id) in enumerate(ids)
            x[j, i] = nodes[node_id][1]
            y[j, i] = nodes[node_id][2]
        end
        x[length(ids)+1, i] = nodes[first(ids)][1]
        y[length(ids)+1, i] = nodes[first(ids)][2]
    end

    seriestype := :path
    linecolor --> :gray
    linewidth --> 0.25
    legend --> false
    aspect_ratio --> :equal
    framestyle --> :box

    return x, y
end
