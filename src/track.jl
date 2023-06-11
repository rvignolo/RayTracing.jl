
@enum DirectionType begin
    Forward
    Backward
end

"""
    Track{BCFwd,BOut,DFwd,DBwd,T<:Real}

Represents a neutron trajectory across the domain, with certain azimuthal angle `ϕ`, entry
and exit points `p` and `q`, length `ℓ` and formed by `segments` coming from its
segmentation.
"""
mutable struct Track{BCFwd,BCBwd,DFwd,DBwd,T<:Real}
    uid::Int
    azim_idx::Int
    track_idx::Int

    p::Point2D{T}
    q::Point2D{T}

    ϕ::T
    ℓ::T

    ABC::SVector{3,T}
    segments::Vector{Segment{T}}

    next_track_fwd::Track
    next_track_bwd::Track

    function Track{BCFwd,BCBwd,DFwd,DBwd}(
        uid::Int,
        azim_idx::Int,
        track_idx::Int,
        p::Point2D{T},
        q::Point2D{T},
        ϕ::T,
        ℓ::T,
        ABC::SVector{3,T},
        segments::Vector{Segment{T}}
    ) where {BCFwd,BCBwd,DFwd,DBwd,T}
        track = new{BCFwd,BCBwd,DFwd,DBwd,T}(
            uid, azim_idx, track_idx, p, q, ϕ, ℓ, ABC, segments
        )
        track.next_track_fwd = track
        track.next_track_bwd = track
        return track
    end
end

universal_id(track::Track) = track.uid
azim_idx(track::Track) = track.azim_idx
track_idx(track::Track) = track.track_idx
bc_fwd(::Track{BCFwd}) where {BCFwd} = BCFwd
bc_bwd(::Track{BCFwd,BCBwd}) where {BCFwd,BCBwd} = BCBwd
dir_next_track_fwd(::Track{BCFwd,BCBwd,DFwd}) where {BCFwd,BCBwd,DFwd} = DFwd
dir_next_track_bwd(::Track{BCFwd,BCBwd,DFwd,DBwd}) where {BCFwd,BCBwd,DFwd,DBwd} = DBwd

function show(io::IO, track::Track)
    @unpack ϕ, p, q, ℓ, segments = track
    # println(io, typeof(t))
    println(io, "  Azimuthal angle: ", round(rad2deg(ϕ), digits=2))
    println(io, "  Entry point: ", p)
    println(io, "  Exit point: ", q)
    println(io, "  Length: ", ℓ)
    println(io, "  # of segments: ", length(segments)) # if it is zero, run segmentize!
    println(io, "  Boundary fwd: ", bc_fwd(track))
    print(io, "  Boundary bwd: ", bc_bwd(track))
    # avoid printing circular references (since it is a cyclic ray tracing...)
    # println(io, "  Next track fwd: ", next_track_fwd)
    # print(io,   "  Next track fwd: ", next_track_bwd)
end

# since we are in active development of the package and there might be unconsidered cases in
# the segmentation of tracks, let's have this maximum number of iterations.
const MAX_ITER = 10_000

function _segmentize_track!(t, track::Track, k::Int, rtol::Real) # t::TrackGenerator
    @unpack mesh, tiny_step = t
    @unpack ϕ, segments = track

    # since we are going to compute them, empty!
    empty!(segments)

    # move a tiny step in ϕ direction to get inside the mesh
    xp = advance_step(track.p, tiny_step, ϕ)

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
                xp = advance_step(xp, tiny_step, ϕ)
                continue
            else
                # or we just finished, so leave
                break
            end
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

        # move a step if the new element is the same as the previous one
        if isequal(prev_element, element)
            xp = advance_step(xp, tiny_step, ϕ)
            continue
        end

        # compute intersections between track and the element
        p, q = intersections(mesh, element, track)

        # we might be in a vertex, so we just continue
        if isapprox(p, q)
            xp = advance_step(xp, tiny_step, ϕ)
            continue
        end

        segment = Segment(p, q, element)
        push!(segments, segment)

        # update new starting point and previous element
        xp = advance_step(q, tiny_step, ϕ)
        prev_element = element

        i += 1 # just for safety
    end

    if !isapprox(track.ℓ, sum(ℓ.(track.segments)); rtol=rtol)
        error("Track with `uid` $(string(track.uid)) has a length that do not match the " *
              "sum of its segments lengths with the provided tolerance `rtol`. Check " *
              "whether this is an actual error or increase `rtol`.")
    end

    return nothing
end
