
@doc raw"""
    general_form(xi::Point2D, xo::Point2D)

Compute the general form equation of a line passing through two points.

This function converts two points into the general form equation `Ax + By + C = 0` that
represents the line passing through both points. The coefficients are normalized for
numerical stability.

## Arguments
- `xi::Point2D`: First point on the line (entry point)
- `xo::Point2D`: Second point on the line (exit point)

## Returns
- `SVector{3,T}`: Coefficients [A, B, C] of the line equation `Ax + By + C = 0`

## Mathematical Details

The general form is computed using the two-point form of a line equation:
- `A = y₁ - y₂` (difference in y-coordinates)
- `B = x₂ - x₁` (difference in x-coordinates, with sign flipped)
- `C = x₁·y₂ - x₂·y₁` (determinant-like term)

The coefficients are then normalized by dividing by `√(A² + B² + C²)` for numerical stability.

## Throws
- `ArgumentError`: If the two points are identical (degenerate line)

## Example
```julia
p1 = Point2D(0.0, 0.0)
p2 = Point2D(1.0, 1.0)
ABC = general_form(p1, p2)  # Returns [A, B, C] for line y = x
```

## Notes

- The function is optimized for performance, avoiding matrix operations
- Normalization prevents numerical issues in subsequent calculations
- Used extensively in ray tracing for track-mesh intersection calculations
"""
function general_form(xi::Point2D, xo::Point2D)
    # Check for degenerate case (identical points)
    if isapprox(xi, xo)
        throw(ArgumentError("Cannot form a line from identical points: $(xi) and $(xo)"))
    end

    # Compute coefficients using direct formula (faster than linear algebra)
    A = xi[2] - xo[2]  # y₁ - y₂
    B = xo[1] - xi[1]  # x₂ - x₁
    C = xi[1] * xo[2] - xo[1] * xi[2]  # x₁·y₂ - x₂·y₁

    # Create coefficient vector and normalize for numerical stability
    # Note: For distinct points, norm cannot be zero mathematically
    ABC = SVector(A, B, C)
    ABC /= norm(ABC)

    return ABC
end

# the function above is faster than using linear algebra:
# [x1 y1; x2 y2] ⋅ [Â, B̂] = [-1, -1] yields to solution for Â and B̂ such that
# Â ⋅ x + B̂ ⋅ y + 1 = 0. Here we make it faster.
# F = transpose(hcat(xi, xo))
# b = SVector(-1, -1)
# x = F \ b
# ABC = vcat(b, 1) # and then we can normalize

"""
    intersections(mesh::Mesh, cell_id::Int, track::Track)

Compute the entry and exit points where a track intersects a mesh element.

This function calculates the intersection points between a ray track and the boundaries of a
specific mesh element. It handles various geometric cases including parallel lines, vertex
intersections, and multiple intersection points.

## Arguments
- `mesh::Mesh`: Computational mesh containing geometry and topology
- `cell_id::Int`: ID of the mesh element to check for intersections
- `track::Track`: Ray track with general form equation `track.ABC`

## Returns
- `Tuple{Point2D, Point2D}`: Entry and exit points `(xi, xo)` of the track through the element

## Algorithm Overview

1. **Element Boundary Traversal**: Iterates through each edge of the mesh element
2. **Line-Line Intersection**: Computes intersection between track and each element edge
3. **Validation**: Filters intersections to ensure they lie on the actual edge segments
4. **Point Selection**: Chooses the appropriate entry/exit points based on intersection count
5. **Ordering**: Orders points according to track direction (azimuthal angle)

## Intersection Cases Handled

### Normal Case (2 intersections)
- Track enters and exits through different edges
- Returns the two intersection points ordered by track direction

### Multiple Intersections (3-4 points)
- Track may intersect at vertices or have complex geometry
- Selects the pair of points with maximum separation distance
- Orders them according to track direction

### Degenerate Cases
- **Parallel lines**: Track parallel to element edge (handled by other edges)
- **Vertex intersection**: Track passes exactly through a vertex
- **Single intersection**: Track grazes the element (returns zero points)
- **No intersections**: Track misses the element (returns zero points)

## Geometric Assumptions

- The track is assumed to intersect the element (caller responsibility)
- Mesh elements are convex polygons (typically triangles or quadrilaterals)
- Track is represented by a straight line segment
- Intersection points are computed with finite precision

## Performance Notes

- Uses pre-allocated `MVector` for intersection points to avoid allocations
- Early termination for parallel lines to avoid unnecessary computations
- Optimized for typical cases (2 intersections) with fallback for edge cases

## Example

```julia
# Compute intersections for a track through element 5
xi, xo = intersections(mesh, 5, track)

# Check if track actually intersects the element
if !isapprox(xi, Point2D(0, 0)) || !isapprox(xo, Point2D(0, 0))
    println("Track enters at $(xi) and exits at $(xo)")
else
    println("No valid intersections found")
end
```

## Notes

- Returns `(Point2D(0,0), Point2D(0,0))` for cases with no valid intersections
- The function is designed for robustness over speed in edge cases
- Used extensively in track segmentation for neutron transport calculations
- Intersection points are ordered according to the track's azimuthal direction

See also: [`intersection`](@ref), [`order_intersection_points`](@ref), [`general_form`](@ref)
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