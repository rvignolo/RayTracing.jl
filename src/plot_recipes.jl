
# TODO: Plots.jl seems slow for these kind of plots (many segments). Use ParaView instead?

@recipe function plot(t::TrackGenerator{T}) where {T}
    @unpack tracks_by_uid, n_total_tracks = t

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
    linewidth --> 0.20
    legend --> false
    framestyle := :none

    return (x, y)
end

@recipe function plot(segments::Vector{Segment{T}}) where {T}

    n = length(segments)
    x = Matrix{T}(undef, 2, n) # TODO: compute Float64
    y = Matrix{T}(undef, 2, n)
    z = Matrix{T}(undef, 2, n)

    for (i, segment) in enumerate(segments)
        x[1, i] = segment.p[1]
        x[2, i] = segment.q[1]
        y[1, i] = segment.p[2]
        y[2, i] = segment.q[2]
        z[1, i] = segment.element
        z[2, i] = segment.element
    end

    seriestype := :path
    linewidth --> 0.25
    legend --> false
    framestyle := :none
    line_z := z

    return (x, y)
end

@recipe function plot(tracks::Vector{Track})

    # m = length(tracks)
    # n = length(segments)
    x = Matrix{Float64}(undef, 2, 0)
    y = Matrix{Float64}(undef, 2, 0)
    z = Matrix{Float64}(undef, 2, 0)

    for track in tracks
        for segment in track.segments
            x = hcat(x, [segment.p[1], segment.q[1]])
            y = hcat(y, [segment.p[2], segment.q[2]])
            z = hcat(z, [segment.element, segment.element])
        end
    end

    seriestype := :path
    linewidth --> 0.25
    legend --> false
    framestyle := :none
    line_z := z

    return (x, y)
end

@recipe function plot(mesh::Mesh)
    @unpack cell_nodes, model = mesh

    grid = get_grid(model)
    nodes = get_node_coordinates(grid)

    # number of cells and number of nodes per cell
    nc = length(cell_nodes)
    nn = all(l -> isequal(l, length(cell_nodes[1])), length.(cell_nodes)) ? length(cell_nodes[1]) : error("error")

    x = Matrix{Float64}(undef, nn + 1, nc) # TODO: compute Float64
    y = Matrix{Float64}(undef, nn + 1, nc)

    for (i, nodes_ids) in enumerate(cell_nodes)

        for (j, node_id) in enumerate(nodes_ids)
            x[j, i] = nodes[node_id][1]
            y[j, i] = nodes[node_id][2]
        end
        # TODO: improve
        x[nn+1, i] = nodes[nodes_ids[1]][1]
        y[nn+1, i] = nodes[nodes_ids[1]][2]
    end

    linecolor --> :black
    seriestype := :path
    linewidth --> 0.2
    legend --> false
    framestyle := :none

    return (x, y)
end