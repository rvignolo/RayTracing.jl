
"""
    TrackGenerator{M<:Mesh,Q<:AzimuthalQuadrature,T<:Real}

This is the main structure of the library, which holds all the information about the ray
tracing.
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

            ABC = general_form(xi, xo)
            segments = Vector{Segment}(undef, 0)

            tracks[i][j] = Track(ϕ, xi, xo, ℓ, ABC, segments)
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

Segmentize tracks, i.e. divide the tracks in segments generated by the intersections between
the cells or elements of the mesh. This function call is intended to be done after calling
`[trace!]`(@ref).
"""
function segmentize!(t::TrackGenerator)
    @unpack tracks_by_uid = t

    !isassigned(tracks_by_uid) && error("Segmentation is intended after tracing. Please, " *
                                        "call `trace!` first!")
    for track in tracks_by_uid
        _segmentize_track!(t, track)
    end
    return t
end