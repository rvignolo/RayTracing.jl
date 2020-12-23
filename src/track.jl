"""
    Track{T<:Real,I<:Real,O<:Real,S}

Represents a neutron trajectory across the domain, with certain azimuthal angle `ϕ`, entry
and exit points `xi` and `xo`, length `ℓ` and formed by `segments` coming from its
segmentation.
"""
struct Track{T<:Real,I<:Point2D,O<:Point2D,S<:AbstractVector,BC<:BoundaryType}
    ϕ::T
    xi::I
    xo::O
    ℓ::T
    ABC::S
    segments::Vector{Segment}

    bi::BC
    bo::BC

    next_track_fwd::Base.RefValue{Track}
    next_track_bwd::Base.RefValue{Track}

    # dir_next_track_fwd::Int
    # dir_next_track_bwd::Int
end

function show(io::IO, track::Track)
    @unpack ϕ, xi, xo, ℓ, segments, bi, bo, next_track_fwd, next_track_bwd = track
    # println(io, typeof(t))
    println(io, "  Azimuthal angle: ", round(rad2deg(ϕ), digits=2))
    println(io, "  Entry point: ", xi)
    println(io, "  Exit point: ", xo)
    println(io, "  Length: ", ℓ)
    println(io, "  # of segments: ", length(segments)) # if it is zero, run segmentize!
    println(io, "  Boundary at entry: ", bi)
    print(io,   "  Boundary at exit: ", bo)
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
    xp = advance_step(track.xi, tiny_step, ϕ)

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
        xi, xo = intersections(mesh, element, track)

        segment = Segment(xi, xo, element)
        push!(segments, segment)

        # update new starting point and previous element
        xp = advance_step(xo, tiny_step, ϕ)
        prev_element = element

        i += 1 # just for safety
    end

    return nothing
end
