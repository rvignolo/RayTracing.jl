
struct TrackGenerator{M,K,B,Q<:AzimuthalQuadrature,T<:Real}
    model::M
    kdtree::K
    bounding_box::B

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

function TrackGenerator(model, n_azim::Int, δ::T; tiny_step::T=1e-8, correct_volumes=false) where {T<:Real}

    n_azim > 0 || throw(DomainError(n_azim, "number of azimuthal angles must be > 0."))
    iszero(rem(n_azim, 4)) || throw(DomainError(n_azim, "number of azimuthal angles must " *
                                                        "be a multiple of 4."))
    δ > 0 || throw(DomainError(δ, "azimuthal spacing must be > 0."))

    grid = get_grid(model)
    kdtree = KDTree(grid)
    bounding_box = BoundingBox(grid)
    Δx, Δy = width(bounding_box), height(bounding_box)

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
        model, kdtree, bounding_box, azimuthal_quadrature, n_tracks_x, n_tracks_y, n_tracks,
        n_total_tracks, tracks, tracks_by_uid, tiny_step, correct_volumes, volumes
    )
end

function trace!(t::TrackGenerator)
    @unpack azimuthal_quadrature, bounding_box = t
    @unpack n_tracks_x, n_tracks_y, n_tracks = t
    @unpack tracks, tracks_by_uid = t
    @unpack n_azim_2, n_azim_4, δs, ϕs, ωₐ = azimuthal_quadrature

    # effective azimuthal spacings, used for some computations but not stored
    T = eltype(δs)
    δx = Vector{T}(undef, n_azim_2)
    δy = Vector{T}(undef, n_azim_2)

    Δx, Δy = width(bounding_box), height(bounding_box)

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
        ωₐ[j] = ωₐ[i]
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
                    xi = Point{2,T}(δx[i] * (n_tracks_x[i] - j + 1 / 2), 0)
                else
                    xi = Point{2,T}(δx[i] * (j - 1 / 2), 0)
                end
            else
                if points_right(azimuthal_quadrature, i)
                    xi = Point{2,T}(0, δy[i] * (j - n_tracks_x[i] - 1 / 2))
                else
                    xi = Point{2,T}(Δx, δy[i] * (j - n_tracks_x[i] - 1 / 2))
                end
            end

            # it can exit at y = Δy regardless if it points to the right or left
            m = tan(ϕ)
            xo = Point{2,T}(xi[1] - (xi[2] - Δy) / m, Δy)

            # keep going if `xo` is not an exit point
            if !(0 <= xo[1] <= Δx)

                # it can exit at x = Δx if it points to the right
                if i <= n_azim_4
                    xo = Point{2,T}(Δx, xi[2] + m * (Δx - xi[1]))

                # or it can exit at x = 0 if it points to the left
                else
                    xo = Point{2,T}(0, xi[2] - m * xi[1])
                end

                if !(0 <= xo[2] <= Δy)
                    throw(DomainError("could not found track exit point."))
                end
            end

            # compute distance and recalibrate coordinates
            ℓ = norm(xi - xo)
            xi = xi + bounding_box.min
            xo = xo + bounding_box.min

            # TODO: mejor calcularlo con LinearAlgebra de alguna forma.
            A = xi[2] - xo[2]
            B = xo[1] - xi[1]
            C = xi[1] * xo[2] - xo[1] * xi[2] # esto pareciera ser un determinante
            ABC = SVector(A, B, C) # por ahi tambien podria ser un Point...
            n = ABC / norm(ABC)

            segments = Vector{Segment}(undef, 0)

            tracks[i][j] = Track(ϕ, xi, xo, ℓ, n, segments)
        end
    end

    uid = 1
    for i in 1:n_azim_2, j in 1:n_tracks[i]
        tracks_by_uid[uid] = tracks[i][j]
        uid += 1
    end

    return t
end

xorigin(n_tracks_x, i, j) = j <= n_tracks_x[i]
yorigin(n_tracks_x, i, j) = !xorigin(n_tracks_x, i, j)

function show(io::IO, t::TrackGenerator{K,B,Q,T}) where {K,B,Q,T}
    @unpack azimuthal_quadrature, bounding_box, n_total_tracks, correct_volumes = t
    @unpack n_azim_2, δs, ϕs = azimuthal_quadrature
    # println(io, typeof(t))
    println(io, "  Number of azimuthal angles in (0, π): ", n_azim_2)
    println(io, "  Azimuthal angles in (0, π): ", round.(rad2deg.(ϕs), digits=2))
    println(io, "  Effective azimuthal spacings: ", round.(δs, digits=3))
    println(io, "  Total tracks: ", n_total_tracks)
    print(io,   "  Correct volumes: ", correct_volumes)
end

@recipe function plot(t::TrackGenerator{K,B,Q,T}; azim_idx=nothing, uid=nothing) where {K,B,Q,T}
    @unpack tracks_by_uid, n_total_tracks = t

    x = Matrix{T}(undef, 2, n_total_tracks)
    y = Matrix{T}(undef, 2, n_total_tracks)

    for j in 1:n_total_tracks
        track = tracks_by_uid[j]
        x[1, j] = track.xi[1]
        x[2, j] = track.xo[1]
        y[1, j] = track.xi[2]
        y[2, j] = track.xo[2]
    end

    # analizar cuales van con --> y con := (una fuerza seguro y la otra no)
    # linecolor   --> :black
    seriestype  :=  :path
    # markershape --> :circle
    linewidth   --> 0.20
    legend      --> false
    border      := :none

    return (x, y)
end


#=
using LazySets
# generate `point`, a vector, and `nodes`, a 2d array or vector of vectors
element = VPolygon(nodes)
println(point in element)
# Also, since you already know that your elements are convex sets and there is no redundance
# in the list of nodes, you can speed up the instantiation by bypassing the convex hull
# algorithm
element = VPolygon(nodes; apply_convex_hull=false)

# oh, an another question: is it possible to compute the intersection between two lines with
# LazySets.jl? Yes, it looks like the Line, Line2D, or LineSegment types will do what you
# need To find the intersection there is a lazy version, the Intersection type, and a
# concrete operation, the intersection function. See
# https://juliareach.github.io/LazySets.jl/dev/lib/binary_functions/#LazySets.intersection-Union{Tuple{N},%20Tuple{Line2D{N,VN}%20where%20VN%3C:AbstractArray{N,1},Line2D{N,VN}%20where%20VN%3C:AbstractArray{N,1}}}%20where%20N%3C:Real
# Disclaimer, I don't know anything about what this package does internally, I'm just a
# happy user for a one-off application where performance is not an issue

=#
