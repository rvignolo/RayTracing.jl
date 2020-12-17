"""
    Track{T<:Real,I<:Real,O<:Real,S}

Represents a neutron trajectory across the domain, with certain azimuthal angle `ϕ`, entry
and exit points `xi` and `xo`, length `ℓ` and formed by `segments` coming from its
segmentation.
"""
struct Track{T<:Real,I<:Point2D,O<:Point2D,S<:AbstractVector}
    ϕ::T
    xi::I
    xo::O
    ℓ::T
    ABC::S
    segments::Vector{Segment}

    # we will eventually need the following fields fields for neutron transport:

    # bi::PhysicalEntity
    # bo::PhysicalEntity

    # next_track_fwd::Track
    # next_track_bwd::Track

    # dir_next_track_fwd::Int
    # dir_next_track_bwd::Int
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

        segment = Segment(xi, xo, norm(xi - xo), element)
        push!(segments, segment)

        # update new starting point and previous element
        xp = advance_step(xo, tiny_step, ϕ)
        prev_element = element

        i += 1 # just for safety
    end

    return nothing
end
