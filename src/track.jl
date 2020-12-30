
@enum DirectionType begin
    Backward
    Forward
end

"""
    Track{UId,AIdx,TIdx,BIn,BOut,DFwd,DBwd,T<:Real}

Represents a neutron trajectory across the domain, with certain azimuthal angle `ϕ`, entry
and exit points `p` and `q`, length `ℓ` and formed by `segments` coming from its
segmentation.
"""
mutable struct Track{UId,AIdx,TIdx,BIn,BOut,DFwd,DBwd,T<:Real}
    p::Point2D{T}
    q::Point2D{T}

    ϕ::T
    ℓ::T

    ABC::SVector{3,T}
    segments::Vector{Segment}

    next_track_fwd::Track
    next_track_bwd::Track

    function Track{UId,AIdx,TIdx,BIn,BOut,DFwd,DBwd}(
        p::Point2D{T}, q::Point2D{T}, ϕ::T, ℓ::T, ABC::SVector{3,T}, segments::Vector{Segment}
    ) where {UId,AIdx,TIdx,BIn,BOut,DFwd,DBwd,T}
        track = new{UId,AIdx,TIdx,BIn,BOut,DFwd,DBwd,T}(p, q, ϕ, ℓ, ABC, segments)
        track.next_track_fwd = track
        track.next_track_bwd = track
        return track
    end
end

universal_id(::Track{UId}) where {UId} = UId
azim_idx(::Track{UId,AIdx}) where {UId,AIdx} = AIdx
track_idx(::Track{UId,AIdx,TIdx}) where {UId,AIdx,TIdx} = TIdx
boundary_in(::Track{UId,AIdx,TIdx,BIn}) where {UId,AIdx,TIdx,BIn} = BIn
boundary_out(::Track{UId,AIdx,TIdx,BIn,BOut}) where {UId,AIdx,TIdx,BIn,BOut} = BOut
dir_next_track_fwd(::Track{UId,AIdx,TIdx,BIn,BOut,DFwd}) where {UId,AIdx,TIdx,BIn,BOut,DFwd} = DFwd
dir_next_track_bwd(::Track{UId,AIdx,TIdx,BIn,BOut,DFwd,DBwd}) where {UId,AIdx,TIdx,BIn,BOut,DFwd,DBwd} = DBwd

function show(io::IO, track::Track)
    @unpack ϕ, p, q, ℓ, segments = track
    # println(io, typeof(t))
    println(io, "  Azimuthal angle: ", round(rad2deg(ϕ), digits=2))
    println(io, "  Entry point: ", p)
    println(io, "  Exit point: ", q)
    println(io, "  Length: ", ℓ)
    println(io, "  # of segments: ", length(segments)) # if it is zero, run segmentize!
    println(io, "  Boundary at entry: ", boundary_in(track))
    print(io,   "  Boundary at exit: ", boundary_out(track))
    # avoid printing circular references (since it is a cyclic ray tracing...)
    # println(io, "  Next track fwd: ", next_track_fwd)
    # print(io,   "  Next track fwd: ", next_track_bwd)
end

# since we are in active development of the package and there might be unconsidered cases in
# the segmentation of tracks, let's have this maximum number of iterations.
const MAX_ITER = 10_000

function _segmentize_track!(t, track::Track, k::Int=5) # t::TrackGenerator
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

        # move a step if the new element is the same as the previous one
        if isequal(prev_element, element)
            xp = advance_step(xp, tiny_step, ϕ)
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
        p, q = intersections(mesh, element, track)

        segment = Segment(p, q, element)
        push!(segments, segment)

        # update new starting point and previous element
        xp = advance_step(q, tiny_step, ϕ)
        prev_element = element

        i += 1 # just for safety
    end

    return nothing
end
