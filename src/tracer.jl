
"""
    TrackGenerator{M<:Mesh,Q<:AzimuthalQuadrature,T<:Real}

This is the main structure of the library, which holds the information of the ray tracing.
"""
struct TrackGenerator{M<:Mesh,Q<:AzimuthalQuadrature,T<:Real}
    mesh::M

    azimuthal_quadrature::Q
    n_tracks_x::Vector{Int}
    n_tracks_y::Vector{Int}
    n_tracks::Vector{Int}
    n_total_tracks::Int

    tracks::Vector{Vector{Track}}
    tracks_by_uid::Vector{Track}

    tiny_step::T
    correct_volumes::Bool
    volumes::Vector{T}
end

function show(io::IO, t::TrackGenerator{M,Q,T}) where {M,Q,T}
    @unpack azimuthal_quadrature, n_total_tracks, correct_volumes = t
    @unpack n_azim_2, δs, ϕs = azimuthal_quadrature
    # println(io, typeof(t))
    println(io, "  Number of azimuthal angles in (0, π): ", n_azim_2)
    println(io, "  Azimuthal angles in (0, π): ", round.(rad2deg.(ϕs), digits=2))
    println(io, "  Effective azimuthal spacings: ", round.(δs, digits=3))
    println(io, "  Total tracks: ", n_total_tracks)
    print(io,   "  Correct volumes: ", correct_volumes)
end

# estas ver donde las mandamos
xorigin(n_tracks_x, i, j) = j <= n_tracks_x[i] # name: origins_in_x
yorigin(n_tracks_x, i, j) = !xorigin(n_tracks_x, i, j)

"""
    TrackGenerator(
        model::UnstructuredDiscreteModel, n_azim::Int, δ::T;
        tiny_step::T=1e-8, correct_volumes=false
    ) where {T<:Real}

Initialize a [`TrackGenerator`](@ref) using an [`UnstructuredDiscreteModel`](@ref) which
holds mesh related information, the number of azimuthal angles `n_azim` and the azimuthal
spacing `δ`. The optional attributes are `tiny_step`, which is used in the track's
segmentation routine, and `correct_volumes`, which allows for cell volume correction since
the ray tracing algorithm yields to approximate volumes.
"""
function TrackGenerator(
    model::UnstructuredDiscreteModel, n_azim::Int, δ::T;
    tiny_step::T=1e-8, correct_volumes=false
) where {T<:Real}

    n_azim > 0 || throw(DomainError(n_azim, "number of azimuthal angles must be > 0."))
    iszero(rem(n_azim, 4)) || throw(DomainError(n_azim, "number of azimuthal angles must " *
                                                        "be a multiple of 4."))
    δ > 0 || throw(DomainError(δ, "azimuthal spacing must be > 0."))

    mesh = Mesh(model)

    Δx, Δy = width(mesh), height(mesh)

    azimuthal_quadrature = AzimuthalQuadrature(n_azim, δ)

    @unpack n_azim_2, n_azim_4 = azimuthal_quadrature

    n_tracks_x = Vector{Int}(undef, n_azim_2)
    n_tracks_y = Vector{Int}(undef, n_azim_2)
    n_tracks   = Vector{Int}(undef, n_azim_2)

    for i in right_dir(azimuthal_quadrature)
        φ = π / n_azim_2 * (i - 1 / 2)

        n_tracks_x[i] = floor(Δx / δ * abs(sin(φ))) + 1
        n_tracks_y[i] = floor(Δy / δ * abs(cos(φ))) + 1
        n_tracks[i] = n_tracks_x[i] + n_tracks_y[i]

        # supplementaries:
        j = suplementary_idx(azimuthal_quadrature, i)
        n_tracks_x[j] = n_tracks_x[i]
        n_tracks_y[j] = n_tracks_y[i]
        n_tracks[j] = n_tracks[i]
    end

    n_total_tracks = sum(n_tracks)

    tracks = Vector{Vector{Track}}(undef, n_azim_2)
    for i in both_dir(azimuthal_quadrature)
        tracks[i] = Vector{Track}(undef, n_tracks[i])
    end

    tracks_by_uid = Vector{Track}(undef, n_total_tracks)

    volumes = Vector{T}(undef, num_cells(model))

    return TrackGenerator(
        mesh, azimuthal_quadrature, n_tracks_x, n_tracks_y, n_tracks, n_total_tracks, tracks,
        tracks_by_uid, tiny_step, correct_volumes, volumes
    )
end

"""
    trace!(t::TrackGenerator)

Computes and fills both the azimuthal quadrature and cyclic tracks around the rectangular
domain using the provided azimuthal angles and spacing when defining the [`TrackGenerator`](@ref)
`t`.
"""
function trace!(t::TrackGenerator{M,Q,T}) where {M,Q,T}
    @unpack mesh, azimuthal_quadrature = t
    @unpack n_tracks_x, n_tracks_y, n_tracks = t
    @unpack tracks, tracks_by_uid = t
    @unpack bbmin, bbmax = mesh
    @unpack n_azim_2, n_azim_4, δs, ϕs, ωₐ = azimuthal_quadrature

    # effective azimuthal spacings, used for some computations but not stored
    δx = Vector{T}(undef, n_azim_2)
    δy = Vector{T}(undef, n_azim_2)

    Δx, Δy = width(mesh), height(mesh)

    for i in right_dir(azimuthal_quadrature)

        # effective azimuthal angle
        ϕ = ϕs[i] = atan((Δy * n_tracks_x[i]) / (Δx * n_tracks_y[i]))

        # effective azimuthal spacings
        δx[i] = Δx / n_tracks_x[i]
        δy[i] = Δy / n_tracks_y[i]
        δs[i] = δx[i] * sin(ϕ)

        # supplementaries:
        j = suplementary_idx(azimuthal_quadrature, i)
        ϕs[j] = π - ϕ
        δx[j] = δx[i]
        δy[j] = δy[i]
        δs[j] = δs[i]
    end

    # once we have computed all the azimuthal angles, compute weights
    init_weights!(azimuthal_quadrature)

    for i in both_dir(azimuthal_quadrature)

        # get azimuthal angle
        ϕ = ϕs[i]

        # for all tracks in a given azimuthal direction
        for j in 1:n_tracks[i]

            if xorigin(n_tracks_x, i, j)
                if points_right(azimuthal_quadrature, i)
                    xi = Point2D(δx[i] * (n_tracks_x[i] - j + 1 / 2), 0)
                else
                    xi = Point2D(δx[i] * (j - 1 / 2), 0)
                end
            else
                if points_right(azimuthal_quadrature, i)
                    xi = Point2D(0, δy[i] * (j - n_tracks_x[i] - 1 / 2))
                else
                    xi = Point2D(Δx, δy[i] * (j - n_tracks_x[i] - 1 / 2))
                end
            end

            # it can exit at y = Δy regardless if it points to the right or left
            m = tan(ϕ)
            xo = Point2D(xi[1] - (xi[2] - Δy) / m, Δy)

            # keep going if `xo` is not an exit point
            if !(0 <= xo[1] <= Δx)

                # it can exit at x = Δx if it points to the right
                if i <= n_azim_4
                    xo = Point2D(Δx, xi[2] + m * (Δx - xi[1]))

                # or it can exit at x = 0 if it points to the left
                else
                    xo = Point2D(0, xi[2] - m * xi[1])
                end

                if !(0 <= xo[2] <= Δy)
                    throw(DomainError("could not found track exit point."))
                end
            end

            # compute distance and recalibrate coordinates
            ℓ = norm(xi - xo)
            xi = xi + bbmin
            xo = xo + bbmin

            n = line_general_equation(xi, xo)
            segments = Vector{Segment}(undef, 0)

            tracks[i][j] = Track(ϕ, xi, xo, ℓ, n, segments)
        end
    end

    # useful pointers for MoC solver
    uid = 1
    for i in 1:n_azim_2, j in 1:n_tracks[i]
        tracks_by_uid[uid] = tracks[i][j]
        uid += 1
    end

    return t
end

"""
    segmentize!(t::TrackGenerator)

Segmentize tracks. This function call is intended to be done after calling `[trace!]`(@ref).
"""
function segmentize!(t::TrackGenerator)
    @unpack tracks_by_uid, tiny_step = t
    for track in tracks_by_uid
        _segmentize_track!(t, track)
    end
    return t
end

const MAX_ITER = 100_000

function _segmentize_track!(t::TrackGenerator, track::Track, k=5)
    @unpack mesh, tiny_step = t
    @unpack ϕ, segments = track

    # since we are going to compute them, empty!
    empty!(segments)

    # move a tiny step in ϕ direction to get inside the mesh
    xp = move_step(track.xi, tiny_step, ϕ)

    i = 0
    element = -1
    prev_element = -1
    while i < MAX_ITER

        # find the element or cell where `xp` lies
        element = find_element(mesh, xp)

        # we might be at the boundary of the mesh
        if inboundary(mesh, xp, tiny_step)
            if isempty(segments)
                # if we just started to segmentize, move a tiny step forward
                xp = move_step(xp, tiny_step, ϕ)
                continue
            else
                # or we just finished, so leave
                break
            end
        end

        # move a step if the new element is the same as the previous one
        if isequal(prev_element, element)
            xp = move_step(xp, tiny_step, ϕ)
            continue
        end

        # if we are inside the domain and `find_element` did not return an element, we might
        # be dealing with a deformed mesh, so we might need to increase the knn search.
        if isequal(element, -1)
            element = find_element(mesh, xp, k)
            if isequal(element, -1)
                error("Try increasing `k`. If the problem persists, raise an issue, this " *
                      "might be a case that hasn't been presented before.")
            end
        end

        # compute intersections between track and the element
        segment_points = intersections(mesh, element, track, tiny_step)

        # TODO: mirar el azim_id del track mejor
        if ϕ < π/2
            if segment_points[1,1] < segment_points[2,1]
                # TODO: son MVector!
                xi = segment_points[1,:]
                xo = segment_points[2,:]
            else
                # TODO: son MVector!
                xi = segment_points[2,:]
                xo = segment_points[1,:]
            end
        else
            if segment_points[1,1] > segment_points[2,1]
                # TODO: son MVector!
                xi = segment_points[1,:]
                xo = segment_points[2,:]
            else
                # TODO: son MVector!
                xi = segment_points[2,:]
                xo = segment_points[1,:]
            end
        end

        xi = convert(Point2D{Float64}, xi)
        xo = convert(Point2D{Float64}, xo)

        segment = Segment(xi, xo, norm(xi - xo), element)
        push!(segments, segment)

        # update new starting point and previous element
        xp = move_step(xo, tiny_step, ϕ)
        prev_element = element

        i += 1 # just for safety
    end

    return nothing
end

@inline move_step(x::Point2D, step::Real, ϕ::Real) = x + step * Point2D(cos(ϕ), sin(ϕ))

function intersections(mesh::Mesh, cell_id::Int, track::Track, tiny_step::Real)
    @unpack model, cell_nodes = mesh

    node_ids = cell_nodes[cell_id]
    node_coordinates = get_node_coordinates(get_grid(model))

    int_points = MMatrix{4,2}(zeros(4,2))
    segment_points = MMatrix{2,2}(zeros(2,2))

    n_int = 0

    parallel_found = 0

    parallel_side = MVector{1,Int64}(0)

    for i in 1:length(node_ids)

        j = i == length(node_ids) ? 1 : i + 1

        x1_id = node_ids[i]
        x2_id = node_ids[j]

        x1 = convert(Point2D{Float64}, node_coordinates[x1_id])
        x2 = convert(Point2D{Float64}, node_coordinates[x2_id])

        ABC = line_general_equation(x1, x2)

        x_int = intersection(track.n, ABC, parallel_side)

        # @show x1
        # @show x2
        # @show x_int

        if parallel_side[] == 1
            parallel_found = 1
            continue
        elseif !point_in_segment(x1, x2, x_int)
            continue
        else
            n_int += 1
            int_points[n_int,:] .= x_int
        end
    end

    # @show n_int
    # @show int_points
    # @show parallel_side[]

    # ahora analizamos en que situacion caimos
    if n_int == 3 || n_int == 4
        ℓ = 0. # TODO: cast zero(T)
        for i in 2:n_int, j in i:n_int

            x1 = int_points[i-1,:]
            x2 = int_points[j,:]

            ℓi = norm(x1 - x2)
            if ℓi > ℓ
                segment_points[1,:] .= x1
                segment_points[2,:] .= x2
                ℓ = ℓi
            end
        end
        return segment_points

    elseif n_int == 2 && parallel_found == 1

        segment_points[1,:] .= int_points[1,:]
        segment_points[2,:] .= int_points[2,:]
        return segment_points

    elseif n_int == 2 && parallel_found == 0

        x1 = int_points[1,:]
        x2 = int_points[2,:]

        if norm(x1 - x2) < tiny_step
            # .... mmmm ... para mi aca tendria que incrementar el tiny step y buscar otro elemento
            @show "estoy aca"
            return compute_intersections(mesh, nodes, track, tiny_step)
        else
            segment_points[1,:] .= int_points[1,:]
            segment_points[2,:] .= int_points[2,:]
            return segment_points
        end
    elseif n_int == 1 || n_int == 0
        return compute_intersections(mesh, nodes, track, tiny_step)
    end
end

function intersection(ABC_1, ABC_2, parallels)

    # calculamos el determinante
    # det = ABC_1[2] * ABC_2[1] - ABC_2[2] * ABC_1[1]

    # si existe (no son paralelas ni coincidentes), computo la interseccion, es decir me
    # fijo que el determinante no sea cero
    if !isapprox(ABC_1[2] * ABC_2[1], ABC_2[2] * ABC_1[1])
        det = ABC_1[2] * ABC_2[1] - ABC_2[2] * ABC_1[1]
        x_int = (ABC_1[3] * ABC_2[2] - ABC_2[3] * ABC_1[2]) / det
        y_int = (ABC_1[1] * ABC_2[3] - ABC_2[1] * ABC_1[3]) / det
        parallels[] = 0
    else
        parallels[] = 1;
    end

    return SVector(x_int, y_int)
end

# TODO: we could just use LinearAlgebra.jl or LazySet.jl
function line_general_equation(xi, xo)
    A = xi[2] - xo[2]
    B = xo[1] - xi[1]
    C = xi[1] * xo[2] - xo[1] * xi[2] # esto pareciera ser un determinante
    ABC = SVector(A, B, C)
    n = ABC / norm(ABC)
    return n
end

# otra forma seria: chequear que el punto cumpla la ecuacion de la recta + inbounds, es
# decir, ya tengo el ABC del segmento, lo uso con el punto x, y luego me fijo que las
# coordenadas de x e y entren dentro del dominio usando inbounds, pero ahi me aparecen
# muchos approx
function point_in_segment(xa, xb, x)
    la = norm(xa - x)
    lb = norm(xb - x)
    ab = norm(xa - xb)
    return isapprox(la + lb, ab)
end
