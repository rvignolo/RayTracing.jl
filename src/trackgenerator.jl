
"""
    TrackGenerator{T<:Real,M<:Mesh,BC<:BoundaryConditions,Q<:AzimuthalQuadrature}

Main structure for ray tracing in unstructured meshes using the Method of Characteristics.

Generates and manages neutron ray trajectories across a 2D domain for transport calculations.
Tracks are organized by azimuthal angle and form closed loops through boundary interactions,
enabling iterative transport sweeps.

## Key Fields
- `mesh`: Computational mesh containing geometry and cell information
- `bcs`: Boundary conditions for domain edges (vacuum, reflective, periodic)
- `azimuthal_quadrature`: Angular discretization with weights and angles
- `tracks`: Tracks organized by azimuthal angle `[azim_idx][track_idx]`
- `tracks_by_uid`: All tracks indexed by unique ID for fast lookup
- `n_total_tracks`: Total number of tracks across all angles
- `tiny_step`: Small step size for numerical stability in intersections
- `volume_correction`: Whether to correct cell volumes from ray tracing
- `volumes`: Cell volumes (corrected if `volume_correction=true`)

## Usage
```julia
# Create track generator
tg = TrackGenerator(model, 16, 0.08, bcs=BoundaryConditions(top=Reflective, ...))

# Generate tracks and segments
trace!(tg)
segmentize!(tg)
```

See also: [`trace!`](@ref), [`segmentize!`](@ref), [`BoundaryConditions`](@ref)
"""
struct TrackGenerator{T<:Real,M<:Mesh,BC<:BoundaryConditions,Q<:AzimuthalQuadrature}
    mesh::M
    bcs::BC

    azimuthal_quadrature::Q
    n_tracks_x::Vector{Int}
    n_tracks_y::Vector{Int}
    n_tracks::Vector{Int}
    n_total_tracks::Int

    tracks::Vector{Vector{Track}}
    tracks_by_uid::Vector{Track}

    tiny_step::T
    volume_correction::Bool
    volumes::Vector{T}
end

function show(io::IO, t::TrackGenerator)
    @unpack azimuthal_quadrature, n_total_tracks, volume_correction = t
    @unpack δs, ϕs = azimuthal_quadrature
    n_azim_2 = nazim2(azimuthal_quadrature)

    # println(io, typeof(t))
    println(io, "  Number of azimuthal angles in (0, π): ", n_azim_2)
    println(io, "  Azimuthal angles in (0, π): ", round.(rad2deg.(ϕs), digits=2))
    println(io, "  Effective azimuthal spacings: ", round.(δs, digits=3))
    println(io, "  Total tracks: ", n_total_tracks)
    print(io, "  Correct volumes: ", volume_correction)
end

origins_in_x(n_tracks_x, i, j) = j ≤ n_tracks_x[i]
origins_in_y(n_tracks_x, i, j) = !origins_in_x(n_tracks_x, i, j)

"""
    TrackGenerator(
        model::UnstructuredDiscreteModel, n_azim::Int, δ::T;
        tiny_step::T=1e-8, volume_correction=false
    ) where {T<:Real}

Initialize a [`TrackGenerator`](@ref) using an [`UnstructuredDiscreteModel`](@ref) which
holds mesh related information, the number of azimuthal angles `n_azim` and the azimuthal
spacing `δ`. The optional attributes are `tiny_step`, which is used in the track's
segmentation routine, and `volume_correction`, which allows for cell volume correction since
the ray tracing algorithm yields to approximate volumes.
"""
function TrackGenerator(
    model::UnstructuredDiscreteModel, n_azim::Int, δ::T;
    bcs::BoundaryConditions=BoundaryConditions(), tiny_step::T=1e-8, volume_correction=false
) where {T<:Real}

    mesh = Mesh(model)
    Δx, Δy = width(mesh), height(mesh)

    azimuthal_quadrature = AzimuthalQuadrature(Val(n_azim), δ)
    n_azim_2 = nazim2(azimuthal_quadrature)
    n_azim_4 = nazim4(azimuthal_quadrature)

    n_tracks_x = Vector{Int}(undef, n_azim_2)
    n_tracks_y = Vector{Int}(undef, n_azim_2)
    n_tracks = Vector{Int}(undef, n_azim_2)

    for i in right_dir(azimuthal_quadrature)
        φ = π / n_azim_2 * (i - 1 / 2)

        n_tracks_x[i] = floor(Δx / δ * abs(sin(φ))) + 1
        n_tracks_y[i] = floor(Δy / δ * abs(cos(φ))) + 1
        n_tracks[i] = n_tracks_x[i] + n_tracks_y[i]

        # suplementary angles:
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
        mesh, bcs, azimuthal_quadrature, n_tracks_x, n_tracks_y, n_tracks, n_total_tracks,
        tracks, tracks_by_uid, tiny_step, volume_correction, volumes
    )
end

"""
    trace!(t::TrackGenerator)

Computes and fills both the azimuthal quadrature and cyclic tracks around the rectangular
domain using the provided azimuthal angles and spacing when defining the [`TrackGenerator`](@ref)
`t`.
"""
function trace!(t::TrackGenerator{T}) where {T}
    @unpack mesh, bcs, azimuthal_quadrature = t
    @unpack n_tracks_x, n_tracks_y, n_tracks = t
    @unpack tracks, tracks_by_uid = t
    @unpack bb_min, bb_max = mesh
    @unpack δs, ϕs, ωₐ = azimuthal_quadrature

    n_azim_2 = nazim2(azimuthal_quadrature)
    n_azim_4 = nazim4(azimuthal_quadrature)

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

        # suplementary angles:
        j = suplementary_idx(azimuthal_quadrature, i)
        ϕs[j] = π - ϕ
        δx[j] = δx[i]
        δy[j] = δy[i]
        δs[j] = δs[i]
    end

    # once we have computed all the azimuthal angles, compute weights
    init_weights!(azimuthal_quadrature)

    # mesh vertices and side segments
    p1 = bb_min
    p2 = Point2D(bb_min[1], bb_max[2])
    p3 = bb_max
    p4 = Point2D(bb_max[1], bb_min[2])
    sides = (top=Segment(p2, p3), bottom=Segment(p4, p1),
        right=Segment(p3, p4), left=Segment(p1, p2))

    uid = 1
    for i in both_dir(azimuthal_quadrature)

        # get azimuthal angle
        ϕ = ϕs[i]

        # for all tracks in a given azimuthal direction `i`
        for j in 1:n_tracks[i]

            if origins_in_x(n_tracks_x, i, j)
                if points_right(azimuthal_quadrature, i)
                    p = Point2D(δx[i] * (n_tracks_x[i] - j + 1 / 2), 0)
                else
                    p = Point2D(δx[i] * (j - 1 / 2), 0)
                end
            else
                if points_right(azimuthal_quadrature, i)
                    p = Point2D(0, δy[i] * (j - n_tracks_x[i] - 1 / 2))
                else
                    p = Point2D(Δx, δy[i] * (j - n_tracks_x[i] - 1 / 2))
                end
            end

            # it can exit at y = Δy regardless if it points to the right or left
            m = tan(ϕ)
            q = Point2D(p[1] - (p[2] - Δy) / m, Δy)

            # keep going if `q` is not an exit point
            if !(0 ≤ q[1] ≤ Δx)

                # it can exit at x = Δx if it points to the right
                if points_right(azimuthal_quadrature, i)
                    q = Point2D(Δx, p[2] + m * (Δx - p[1]))

                else
                    # or it can exit at x = 0 if it points to the left
                    q = Point2D(0, p[2] - m * p[1])
                end

                if !(0 ≤ q[2] ≤ Δy)
                    throw(DomainError("could not found track exit point."))
                end
            end

            # recalibrate coordinates and compute distance
            p += bb_min
            q += bb_min
            ℓ = norm(p - q)

            ABC = general_form(p, q)
            segments = Vector{Segment{T}}(undef, 0)

            BCFwd = boundary_condition(q, sides, bcs)
            BCBwd = boundary_condition(p, sides, bcs)

            # alternative method
            if points_right(azimuthal_quadrature, i)
                BCFwd1 = j ≤ n_tracks_y[i] ? bcs.right : bcs.top
                BCBwd1 = j ≤ n_tracks_x[i] ? bcs.bottom : bcs.left
            else
                BCFwd1 = j ≤ n_tracks_y[i] ? bcs.left : bcs.top
                BCBwd1 = j ≤ n_tracks_x[i] ? bcs.bottom : bcs.right
            end

            if !isequal(BCFwd, BCFwd1) || !isequal(BCBwd, BCBwd1)
                error("Boundaries do not match!")
            end

            if j ≤ n_tracks_y[i]
                DirNextTrackFwd = Forward
            else
                if BCFwd == Periodic
                    DirNextTrackFwd = Forward
                elseif BCFwd == Vacuum || BCFwd == Reflective
                    DirNextTrackFwd = Backward
                end
            end

            if j ≤ n_tracks_x[i]
                if BCBwd == Periodic
                    DirNextTrackBwd = Backward
                elseif BCBwd == Vacuum || BCBwd == Reflective
                    DirNextTrackBwd = Forward
                end
            else
                DirNextTrackBwd = Backward
            end

            track = Track{BCFwd,BCBwd,DirNextTrackFwd,DirNextTrackBwd}(
                uid, i, j, p, q, ϕ, ℓ, ABC, segments
            )
            tracks_by_uid[uid] = track
            tracks[i][j] = track
            uid += 1
        end
    end

    # once we have all tracks and boundary conditions, set next tracks
    next_tracks(t)

    return t
end

function next_tracks(t::TrackGenerator)
    @unpack tracks_by_uid = t

    # search for the next track in forward and backward direction
    for track in tracks_by_uid
        next_track_fwd(t, track)
        next_track_bwd(t, track)
    end

    return nothing
end

function next_track_fwd(t::TrackGenerator, track::Track)
    @unpack azimuthal_quadrature, tracks = t
    @unpack n_tracks_x, n_tracks_y, n_tracks = t
    @unpack azim_idx, track_idx = track

    i, j, k = azim_idx, track_idx, suplementary_idx(azimuthal_quadrature, azim_idx)

    BCFwd = bc_fwd(track)

    # these are the tracks that arrive to the y-axis
    if j ≤ n_tracks_y[i]
        if BCFwd == Periodic
            track.next_track_fwd = tracks[i][j+n_tracks_x[i]]
        elseif BCFwd == Vacuum || BCFwd == Reflective
            track.next_track_fwd = tracks[k][j+n_tracks_x[i]]
        end
    else
        # these are the tracks that arrive to the top (superior x-axis)
        if BCFwd == Periodic
            track.next_track_fwd = tracks[i][j-n_tracks_y[i]]
        elseif BCFwd == Vacuum || BCFwd == Reflective
            track.next_track_fwd = tracks[k][n_tracks[i]+n_tracks_y[i]-j+1]
        end
    end

    return nothing
end

function next_track_bwd(t::TrackGenerator, track::Track)
    @unpack azimuthal_quadrature, tracks = t
    @unpack n_tracks_x, n_tracks_y, n_tracks = t
    @unpack azim_idx, track_idx = track

    i, j, k = azim_idx, track_idx, suplementary_idx(azimuthal_quadrature, azim_idx)

    BCBwd = bc_bwd(track)

    # these are the tracks that arrive to the bottom (inferior x-axis)
    if j ≤ n_tracks_x[i]
        if BCBwd == Periodic
            track.next_track_bwd = tracks[i][j+n_tracks_y[i]]
        elseif BCBwd == Vacuum || BCBwd == Reflective
            track.next_track_bwd = tracks[k][n_tracks_x[i]-j+1]
        end
        # these are the tracks that arrive to the y-axis
    else
        if BCBwd == Periodic
            track.next_track_bwd = tracks[i][j-n_tracks_x[i]]
        elseif BCBwd == Vacuum || BCBwd == Reflective
            track.next_track_bwd = tracks[k][j-n_tracks_x[i]]
        end
    end

    return nothing
end

"""
    segmentize!(t::TrackGenerator)

Segmentize tracks, i.e. divide the tracks in segments generated by the intersections between
the cells or elements of the mesh. This function call is intended to be done after calling
`[trace!]`(@ref).
"""
function segmentize!(t::TrackGenerator{T}; k::Int=5, rtol::Real=Base.rtoldefault(T)) where {T}
    @unpack tracks_by_uid = t

    !isassigned(tracks_by_uid, 1) && error("Segmentation is intended after tracing. Please, " *
                                           "call `trace!` first!")
    for track in tracks_by_uid
        _segmentize_track!(t, track, k, rtol)
    end

    fill_volumes(t, 1)

    return t
end

function fill_volumes(t::TrackGenerator{T}, i) where {T}
    @unpack tracks_by_uid, azimuthal_quadrature, volumes = t
    @unpack δs = azimuthal_quadrature
    n_azim_2 = nazim2(azimuthal_quadrature)

    fill!(volumes, zero(T))

    for track in tracks_by_uid
        a = track.azim_idx
        for segment in track.segments
            i = segment.element
            volumes[i] += δs[a] * segment.ℓ
        end
    end

    volumes ./= n_azim_2

    #### TODO: correct volumes by changing segment lengths #####
    @unpack mesh = t
    @unpack model, cell_nodes = mesh
    ncells = num_cells(model)
    volumes2 = Vector{T}(undef, ncells)
    fill!(volumes2, zero(T))
    for i in 1:ncells
        node_ids = cell_nodes[i]
        volumes2[i] = element_volume(mesh, node_ids)
    end

    return nothing
end

function element_volume(mesh, node_ids)
    @unpack model = mesh
    node_coordinates = get_node_coordinates(get_grid(model))

    x1 = node_coordinates[node_ids[1]]
    x2 = node_coordinates[node_ids[2]]
    x3 = node_coordinates[node_ids[3]]

    return 1 / 2 * abs((x2 - x1) × (x3 - x1))
end
