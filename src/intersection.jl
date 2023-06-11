
@doc raw"""
    general_form(xi::Point2D, xo::Point2D)

Returns the general form equation of a line that passes through points `xi` and `xo` as a
[`SVector`](@ref) which holds A, B and C such that:
```math
A \cdot x + B \cdot y +  C = 0.
```
"""
function general_form(xi::Point2D, xo::Point2D)
    A = xi[2] - xo[2]
    B = xo[1] - xi[1]
    C = xi[1] * xo[2] - xo[1] * xi[2]
    ABC = SVector(A, B, C)
    ABC /= norm(ABC) # for numerical reasons (?)
    return ABC
end
# the function above is faster than using linear algebra:
# [x1 y1; x2 y2] ⋅ [Â, B̂] = [-1, -1] yields to solution for Â and B̂ such that
# Â ⋅ x + B̂ ⋅ y + 1 = 0. Here we make it faster.
# F = transpose(hcat(xi, xo))
# b = SVector(-1, -1)
# x = F \ b
# ABC = vcat(b, 1) # and then we can normalize

# We could use LazySets.jl for some part of this particular function but it seems slow!
"""
    intersections(mesh::Mesh, cell_id::Int, track::Track)

Computes the entry and exit points of a given `track` that is known to cross element with id
`cell_id`.
"""
function intersections(
    mesh::Mesh, cell_id::Int32, track::Track{BIn,BOut,DFwd,DBwd,T}
) where {BIn,BOut,DFwd,DBwd,T}
    @unpack model, cell_nodes = mesh

    node_coordinates = get_node_coordinates(get_grid(model))
    node_ids = cell_nodes[cell_id]

    int_points = MVector{4,Point2D{T}}(
        Point2D(0, 0), Point2D(0, 0), Point2D(0, 0), Point2D(0, 0)
    )

    n_int = 0
    parallel_found = false
    for i in eachindex(node_ids)

        j = i == lastindex(node_ids) ? 1 : i + 1

        # get node coordinates and cast them to Point2D
        p1 = convert(Point2D{T}, node_coordinates[node_ids[i]])
        p2 = convert(Point2D{T}, node_coordinates[node_ids[j]])

        # compute general form equation for the selected element face
        ABC = general_form(p1, p2)

        # compute intersections between track and element face
        are_parallel, x_int = intersection(track.ABC, ABC)

        # if the track is parallel to the element face, it can be on top of it or do not
        # cross at all. Either way, we just ignore this case because we can compute the
        # intersections using the other faces.
        if are_parallel
            parallel_found = true
            continue

        elseif !point_in_segment(p1, p2, x_int)
            # if the intersection is outside the face, avoid it
            continue

        else
            # this is a valid intersection, store it
            n_int += 1
            int_points[n_int] = x_int
        end
    end

    # now let's check all the possible situations and handle extreme cases
    if n_int in (3, 4)
        ℓ = zero(T)
        # get the points that have the maximun distance between
        for i in 2:n_int, j in i:n_int
            x1 = int_points[i-1]
            x2 = int_points[j]
            ℓi = norm(x1 - x2)
            if ℓi > ℓ
                x_int1 = x1
                x_int2 = x2
                ℓ = ℓi
            end
        end
        return order_intersection_points(track, x_int1, x_int2)

    elseif n_int == 2 && parallel_found

        x_int1 = int_points[1]
        x_int2 = int_points[2]
        return order_intersection_points(track, x_int1, x_int2)

    elseif n_int == 2 && !parallel_found

        x_int1 = int_points[1]
        x_int2 = int_points[2]

        if isapprox(x_int1, x_int2)
            # do nothing and move a tiny step forward in the parent function, we are on a
            # vertex
            return x_int1, x_int2
        else
            return order_intersection_points(track, x_int1, x_int2)
        end
    elseif iszero(n_int) || isone(n_int)
        # the parent function will move a tiny step further
        return Point2D{T}(0, 0), Point2D{T}(0, 0)
        # error("This is an unexpected case. Please, submit an issue.")
    end
end

"""
    intersection(ABC1, ABC2)

Computes the intersection between two lines with general form equations given by `ABC1` and
`ABC2`. It also indicates if the lines are parallel.
"""
function intersection(ABC1::AbstractVector, ABC2::AbstractVector)
    a = ABC1[2] * ABC2[1]
    b = ABC2[2] * ABC1[1]
    x = y = zero(a)
    are_parallel = isapprox(a, b)
    if !are_parallel
        det = a - b
        x = (ABC1[3] * ABC2[2] - ABC2[3] * ABC1[2]) / det
        y = (ABC1[1] * ABC2[3] - ABC2[1] * ABC1[3]) / det
    end
    return are_parallel, Point2D(x, y)
end

# this one returns a matrix even if F is not invertible (has the same speed as the above)
# function intersection(ABC1, ABC2)
#     F = @SMatrix [ABC1[1] ABC1[2]
#                   ABC2[1] ABC2[2]]
#     b = @SVector [ABC1[3], ABC2[3]]

#     x = F \ (-b)

#     return Point2D(x[1], x[2])
# end

function order_intersection_points(track::Track, x1::Point2D, x2::Point2D)
    @unpack ϕ = track
    if isless(ϕ, π / 2)
        xi, xo = first(x1) < first(x2) ? (x1, x2) : (x2, x1)
    else
        xi, xo = first(x1) > first(x2) ? (x1, x2) : (x2, x1)
    end
    return xi, xo
end