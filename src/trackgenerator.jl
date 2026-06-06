
"""
    TrackGenerator{T<:Real,M<:Mesh,BC<:BoundaryConditions,Q<:AzimuthalQuadrature}

Main structure for ray tracing in unstructured meshes using the Method of Characteristics.

Generates and manages neutron ray trajectories across a 2D domain for transport calculations. Tracks
are organized by azimuthal angle and form closed loops through boundary interactions, enabling
iterative transport sweeps.

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
# Create a track generator.
tg = TrackGenerator(model, 16, 0.08, bcs=BoundaryConditions(top=Reflective, ...))

# Generate tracks and segments.
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

    tracks::Vector{Vector{Track{T}}}
    tracks_by_uid::Vector{Track{T}}

    tiny_step::T
    volume_correction::Bool
    volumes::Vector{T}
end

function show(io::IO, t::TrackGenerator)
    @unpack azimuthal_quadrature, n_total_tracks, volume_correction = t
    @unpack δs, ϕs = azimuthal_quadrature
    n_azim_2 = n_azim_half(azimuthal_quadrature)

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

Initialize a [`TrackGenerator`](@ref) using an [`UnstructuredDiscreteModel`](@ref), the number
of azimuthal angles `n_azim`, and the azimuthal spacing `δ`. The optional attributes are
`tiny_step`, which is used in track segmentation, and `volume_correction`, which corrects cell
volumes because ray tracing produces approximate volumes.
"""
function TrackGenerator(
    model::UnstructuredDiscreteModel, n_azim::Int, δ::T;
    bcs::BoundaryConditions=BoundaryConditions(), tiny_step::T=1e-8, volume_correction=false
) where {T<:Real}

    mesh = Mesh(model)
    Δx, Δy = width(mesh), height(mesh)

    azimuthal_quadrature = AzimuthalQuadrature(Val(n_azim), δ)
    n_azim_2 = n_azim_half(azimuthal_quadrature)

    n_tracks_x = Vector{Int}(undef, n_azim_2)
    n_tracks_y = Vector{Int}(undef, n_azim_2)
    n_tracks = Vector{Int}(undef, n_azim_2)

    for i in azimuthal_quadrant_1(azimuthal_quadrature)
        φ = π / n_azim_2 * (i - 1 / 2)

        n_tracks_x[i] = floor(Δx / δ * abs(sin(φ))) + 1
        n_tracks_y[i] = floor(Δy / δ * abs(cos(φ))) + 1
        n_tracks[i] = n_tracks_x[i] + n_tracks_y[i]

        # Supplementary angles:
        j = supplementary_azimuthal_idx(azimuthal_quadrature, i)
        n_tracks_x[j] = n_tracks_x[i]
        n_tracks_y[j] = n_tracks_y[i]
        n_tracks[j] = n_tracks[i]
    end

    n_total_tracks = sum(n_tracks)

    tracks = Vector{Vector{Track{T}}}(undef, n_azim_2)
    for i in azimuthal_half_plane(azimuthal_quadrature)
        tracks[i] = Vector{Track{T}}(undef, n_tracks[i])
    end

    tracks_by_uid = Vector{Track{T}}(undef, n_total_tracks)

    volumes = Vector{T}(undef, num_cells(model))

    return TrackGenerator(
        mesh, bcs, azimuthal_quadrature, n_tracks_x, n_tracks_y, n_tracks, n_total_tracks,
        tracks, tracks_by_uid, tiny_step, volume_correction, volumes
    )
end

"""
    trace!(t::TrackGenerator)

Computes and fills both the azimuthal quadrature and cyclic tracks around the rectangular domain
using the provided azimuthal angles and spacing when defining the [`TrackGenerator`](@ref) `t`.
"""
function trace!(t::TrackGenerator{T}) where {T}
    @unpack mesh, bcs, azimuthal_quadrature = t
    @unpack n_tracks_x, n_tracks_y, n_tracks = t
    @unpack tracks, tracks_by_uid = t
    @unpack bb_min, bb_max = mesh
    @unpack δs, ϕs = azimuthal_quadrature

    n_azim_2 = n_azim_half(azimuthal_quadrature)

    # Effective azimuthal spacings used for intermediate computations.
    δx = Vector{T}(undef, n_azim_2)
    δy = Vector{T}(undef, n_azim_2)

    Δx, Δy = width(mesh), height(mesh)

    for i in azimuthal_quadrant_1(azimuthal_quadrature)

        # Effective azimuthal angle.
        ϕ = ϕs[i] = atan((Δy * n_tracks_x[i]) / (Δx * n_tracks_y[i]))

        # Effective azimuthal spacings.
        δx[i] = Δx / n_tracks_x[i]
        δy[i] = Δy / n_tracks_y[i]
        δs[i] = δx[i] * sin(ϕ)

        # Supplementary angles:
        j = supplementary_azimuthal_idx(azimuthal_quadrature, i)
        ϕs[j] = π - ϕ
        δx[j] = δx[i]
        δy[j] = δy[i]
        δs[j] = δs[i]
    end

    # Compute weights after all azimuthal angles are available.
    init_weights!(azimuthal_quadrature)

    # Mesh vertices and boundary segments.
    p1 = bb_min
    p2 = Point2D(bb_min[1], bb_max[2])
    p3 = bb_max
    p4 = Point2D(bb_max[1], bb_min[2])
    boundary = Boundary(Segment(p2, p3), Segment(p4, p1), Segment(p3, p4), Segment(p1, p2))

    uid = 1
    for i in azimuthal_half_plane(azimuthal_quadrature)

        # Get the azimuthal angle.
        ϕ = ϕs[i]

        # Iterate over all tracks in azimuthal direction `i`.
        for j in 1:n_tracks[i]

            if origins_in_x(n_tracks_x, i, j)
                if is_rightward_direction(azimuthal_quadrature, i)
                    p = Point2D(δx[i] * (n_tracks_x[i] - j + 1 / 2), 0)
                else
                    p = Point2D(δx[i] * (j - 1 / 2), 0)
                end
            else
                if is_rightward_direction(azimuthal_quadrature, i)
                    p = Point2D(0, δy[i] * (j - n_tracks_x[i] - 1 / 2))
                else
                    p = Point2D(Δx, δy[i] * (j - n_tracks_x[i] - 1 / 2))
                end
            end

            # The track can exit at y = Δy whether it points right or left.
            m = tan(ϕ)
            q = Point2D(p[1] - (p[2] - Δy) / m, Δy)

            # Continue searching if `q` is not a valid exit point.
            if !(0 ≤ q[1] ≤ Δx)

                # A rightward track can exit at x = Δx.
                if is_rightward_direction(azimuthal_quadrature, i)
                    q = Point2D(Δx, p[2] + m * (Δx - p[1]))

                else
                    # A leftward track can exit at x = 0.
                    q = Point2D(0, p[2] - m * p[1])
                end

                if !(0 ≤ q[2] ≤ Δy)
                    throw(DomainError("could not found track exit point."))
                end
            end

            # Shift coordinates into the mesh frame and compute distance.
            p += bb_min
            q += bb_min
            ℓ = norm(p - q)

            ABC = general_form(p, q)
            segments = Vector{Segment{T}}(undef, 0)

            BCFwd = get_boundary_condition_at(q, boundary, bcs)
            BCBwd = get_boundary_condition_at(p, boundary, bcs)

            # Cross-check boundary conditions using track indices.
            if is_rightward_direction(azimuthal_quadrature, i)
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
                if is_periodic(BCFwd)
                    DirNextTrackFwd = Forward
                elseif is_vacuum(BCFwd) || is_reflective(BCFwd)
                    DirNextTrackFwd = Backward
                end
            end

            if j ≤ n_tracks_x[i]
                if is_periodic(BCBwd)
                    DirNextTrackBwd = Backward
                elseif is_vacuum(BCBwd) || is_reflective(BCBwd)
                    DirNextTrackBwd = Forward
                end
            else
                DirNextTrackBwd = Backward
            end

            track = Track(
                uid, i, j, BCFwd, BCBwd, DirNextTrackFwd, DirNextTrackBwd,
                p, q, ϕ, ℓ, ABC, segments
            )
            tracks_by_uid[uid] = track
            tracks[i][j] = track
            uid += 1
        end
    end

    # Set next-track links after all tracks and boundary conditions are known.
    next_tracks(t)

    return t
end

function next_tracks(t::TrackGenerator)
    @unpack tracks_by_uid = t

    # Search for next tracks in the forward and backward directions.
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

    i, j, k = azim_idx, track_idx, supplementary_azimuthal_idx(azimuthal_quadrature, azim_idx)

    BCFwd = bc_fwd(track)

    # Tracks that arrive at the y-axis.
    if j ≤ n_tracks_y[i]
        if BCFwd == Periodic
            track.next_track_fwd = tracks[i][j+n_tracks_x[i]]
        elseif is_vacuum(BCFwd) || is_reflective(BCFwd)
            track.next_track_fwd = tracks[k][j+n_tracks_x[i]]
        end
    else
        # Tracks that arrive at the top edge.
        if BCFwd == Periodic
            track.next_track_fwd = tracks[i][j-n_tracks_y[i]]
        elseif is_vacuum(BCFwd) || is_reflective(BCFwd)
            track.next_track_fwd = tracks[k][n_tracks[i]+n_tracks_y[i]-j+1]
        end
    end

    return nothing
end

function next_track_bwd(t::TrackGenerator, track::Track)
    @unpack azimuthal_quadrature, tracks = t
    @unpack n_tracks_x, n_tracks_y, n_tracks = t
    @unpack azim_idx, track_idx = track

    i, j, k = azim_idx, track_idx, supplementary_azimuthal_idx(azimuthal_quadrature, azim_idx)

    BCBwd = bc_bwd(track)

    # Tracks that arrive at the bottom edge.
    if j ≤ n_tracks_x[i]
        if is_periodic(BCBwd)
            track.next_track_bwd = tracks[i][j+n_tracks_y[i]]
        elseif is_vacuum(BCBwd) || is_reflective(BCBwd)
            track.next_track_bwd = tracks[k][n_tracks_x[i]-j+1]
        end
        # Tracks that arrive at the y-axis.
    else
        if is_periodic(BCBwd)
            track.next_track_bwd = tracks[i][j-n_tracks_x[i]]
        elseif is_vacuum(BCBwd) || is_reflective(BCBwd)
            track.next_track_bwd = tracks[k][j-n_tracks_x[i]]
        end
    end

    return nothing
end

"""
    segmentize!(t::TrackGenerator; parallel=:auto)

Split tracks into segments generated by intersections with mesh cells/elements. Call this
after [`trace!`](@ref). The `parallel` keyword accepts `:auto`, `true`, or `false`.
With `parallel=:auto`, segmentation uses `Threads.@threads` when Julia has multiple threads.
Volume accumulation is still performed after all track segments are computed.
"""
function segmentize!(
    t::TrackGenerator{T}; k::Int=5, rtol::Real=Base.rtoldefault(T), parallel=:auto
) where {T}
    @unpack tracks_by_uid = t

    !isassigned(tracks_by_uid, 1) && error("Segmentation is intended after tracing. Please, " *
                                           "call `trace!` first!")
    if should_segmentize_parallel(parallel)
        segmentize_tracks_threaded!(t, k, rtol)
    else
        segmentize_tracks_serial!(t, k, rtol)
    end

    fill_volumes(t)

    return t
end

function should_segmentize_parallel(parallel)
    parallel === false && return false
    parallel === true && return Threads.nthreads() > 1
    parallel === :auto && return Threads.nthreads() > 1
    throw(ArgumentError("parallel must be :auto, true, or false; got $parallel"))
end

function segmentize_tracks_serial!(t::TrackGenerator, k::Int, rtol::Real)
    @unpack tracks_by_uid = t
    for track in tracks_by_uid
        _segmentize_track!(t, track, k, rtol)
    end
    return nothing
end

function segmentize_tracks_threaded!(t::TrackGenerator, k::Int, rtol::Real)
    @unpack tracks_by_uid = t
    Threads.@threads for i in eachindex(tracks_by_uid)
        _segmentize_track!(t, tracks_by_uid[i], k, rtol)
    end
    return nothing
end

function fill_volumes(t::TrackGenerator{T}) where {T}
    @unpack mesh, tracks_by_uid, azimuthal_quadrature, volumes, volume_correction = t
    @unpack δs = azimuthal_quadrature
    n_azim_2 = n_azim_half(azimuthal_quadrature)

    fill!(volumes, zero(T))

    for track in tracks_by_uid
        accumulate_track_volume!(volumes, δs, track)
    end

    volumes ./= n_azim_2

    if volume_correction
        for i in eachindex(volumes)
            volumes[i] = element_volume(mesh, i)
        end
    end

    return nothing
end

function accumulate_track_volume!(volumes, δs, track::Track)
    δ = δs[track.azim_idx]
    for segment in track.segments
        volumes[segment.element] += δ * segment.ℓ
    end
    return nothing
end

function element_volume(mesh, node_ids)
    @unpack node_coordinates = mesh
    ordered_ids = ordered_node_ids(mesh, node_ids)
    return _element_volume_from_ordered_nodes(node_coordinates, ordered_ids)
end

function element_volume(mesh::Mesh, cell_id::Integer)
    @unpack node_coordinates, ordered_cell_nodes = mesh
    return _element_volume_from_ordered_nodes(node_coordinates, ordered_cell_nodes[cell_id])
end

function _element_volume_from_ordered_nodes(node_coordinates, ordered_ids)
    area = zero(eltype(first(node_coordinates)))
    for i in eachindex(ordered_ids)
        j = i == lastindex(ordered_ids) ? firstindex(ordered_ids) : i + 1
        p = node_coordinates[ordered_ids[i]]
        q = node_coordinates[ordered_ids[j]]
        area += p[1] * q[2] - q[1] * p[2]
    end

    return abs(area) / 2
end
