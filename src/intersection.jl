
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

## Notes

- Returns `(Point2D(0,0), Point2D(0,0))` for cases with no valid intersections
- The function is designed for robustness over speed in edge cases
- Used extensively in track segmentation for neutron transport calculations
- Intersection points are ordered according to the track's azimuthal direction

See also: [`intersection`](@ref), [`order_intersection_points`](@ref), [`general_form`](@ref)
"""
function intersections(
    mesh::Mesh, cell_id::Int32, track::Track{BCFwd,BCBwd,DFwd,DBwd,T}
) where {BCFwd,BCBwd,DFwd,DBwd,T}

    # Validate inputs
    cell_id > 0 || throw(ArgumentError("Cell ID must be positive, got: $cell_id"))
    cell_id <= length(mesh.cell_nodes) || throw(ArgumentError("Cell ID out of range: $cell_id"))

    @unpack model, cell_nodes = mesh
    node_coordinates = get_node_coordinates(get_grid(model))
    cell_node_ids = cell_nodes[cell_id]

    # Validate node_ids
    length(cell_node_ids) >= 3 || throw(ArgumentError("Element must have at least 3 nodes"))

    intersection_points = MVector{4,Point2D{T}}(
        Point2D(0, 0), Point2D(0, 0), Point2D(0, 0), Point2D(0, 0)
    )

    num_intersections = 0
    has_parallel_edge = false

    for edge_idx in eachindex(cell_node_ids)

        next_edge_idx = edge_idx == lastindex(cell_node_ids) ? 1 : edge_idx + 1

        # get node coordinates and cast them to Point2D
        p1 = convert(Point2D{T}, node_coordinates[cell_node_ids[edge_idx]])
        p2 = convert(Point2D{T}, node_coordinates[cell_node_ids[next_edge_idx]])

        # compute general form equation for the selected element face
        ABC = general_form(p1, p2)

        # compute intersections between track and element face
        are_parallel, x_int = intersection(track.ABC, ABC)

        # if the track is parallel to the element face, it can be on top of it or do not
        # cross at all. Either way, we just ignore this case because we can compute the
        # intersections using the other faces.
        if are_parallel
            has_parallel_edge = true
            continue

        elseif !point_in_segment(p1, p2, x_int)
            # if the intersection is outside the face, avoid it
            continue

        else
            # this is a valid intersection, store it
            num_intersections += 1
            intersection_points[num_intersections] = x_int
        end
    end

    # now let's check all the possible situations and handle extreme cases
    if num_intersections in (3, 4)
        ℓ = zero(T)
        # get the points that have the maximun distance between
        for i in 2:num_intersections, j in i:num_intersections
            x1 = intersection_points[i-1]
            x2 = intersection_points[j]
            ℓi = norm(x1 - x2)
            if ℓi > ℓ
                x_int1 = x1
                x_int2 = x2
                ℓ = ℓi
            end
        end
        return order_intersection_points(track, x_int1, x_int2)

    elseif num_intersections == 2 && has_parallel_edge

        x_int1 = intersection_points[1]
        x_int2 = intersection_points[2]
        return order_intersection_points(track, x_int1, x_int2)

    elseif num_intersections == 2 && !has_parallel_edge

        x_int1 = intersection_points[1]
        x_int2 = intersection_points[2]

        if isapprox(x_int1, x_int2)
            # do nothing and move a tiny step forward in the parent function, we are on a
            # vertex
            return x_int1, x_int2
        else
            return order_intersection_points(track, x_int1, x_int2)
        end
    elseif iszero(num_intersections) || isone(num_intersections)
        # the parent function will move a tiny step further
        return Point2D{T}(0, 0), Point2D{T}(0, 0)
        # error("This is an unexpected case. Please, submit an issue.")
    end
end

"""
    intersection(ABC1::AbstractVector, ABC2::AbstractVector)

Compute the intersection point between two lines given in general form.

This function finds the intersection point of two lines represented by their general form
equations `A₁x + B₁y + C₁ = 0` and `A₂x + B₂y + C₂ = 0`. It also detects when the lines
are parallel (no intersection).

## Arguments
- `ABC1::AbstractVector`: Coefficients [A₁, B₁, C₁] of the first line equation
- `ABC2::AbstractVector`: Coefficients [A₂, B₂, C₂] of the second line equation

## Returns
- `Tuple{Bool, Point2D}`: `(are_parallel, intersection_point)`
  - `are_parallel`: `true` if lines are parallel, `false` otherwise
  - `intersection_point`: Point of intersection (meaningful only if `are_parallel = false`)

## Mathematical Details

The intersection is computed by solving the system of linear equations:
```
A₁x + B₁y = -C₁
A₂x + B₂y = -C₂
```

Using Cramer's rule:
- **Determinant**: `det = A₁B₂ - A₂B₁`
- **Parallel lines**: `det ≈ 0` (lines are parallel or coincident)
- **Intersection point**:
  - `x = (B₁C₂ - B₂C₁) / det`
  - `y = (A₂C₁ - A₁C₂) / det`

## Example

```julia
# Line 1: x + y - 1 = 0
ABC1 = [1.0, 1.0, -1.0]

# Line 2: x - y - 1 = 0
ABC2 = [1.0, -1.0, -1.0]

are_parallel, point = intersection(ABC1, ABC2)
# Returns: (false, Point2D(1.0, 0.0))
```

## Notes

- Uses `isapprox` for parallel detection to handle numerical precision issues
- Returns `Point2D(0, 0)` when lines are parallel (convention for no intersection)
- The function is optimized for performance in ray tracing applications
- Coefficients should be normalized for best numerical stability

See also: [`general_form`](@ref), [`intersections`](@ref)
"""
function intersection(ABC1::AbstractVector, ABC2::AbstractVector)
    # Extract coefficients for clarity
    A1, B1, C1 = ABC1[1], ABC1[2], ABC1[3]
    A2, B2, C2 = ABC2[1], ABC2[2], ABC2[3]

    # Compute determinant for parallel detection
    determinant = A1 * B2 - A2 * B1

    # Check if lines are parallel (determinant ≈ 0)
    z = zero(determinant)
    are_parallel = isapprox(determinant, z)

    if are_parallel
        # Lines are parallel or coincident - no unique intersection
        return true, Point2D(z, z)
    else
        # Compute intersection point using Cramer's rule
        x = (B1 * C2 - B2 * C1) / determinant
        y = (A2 * C1 - A1 * C2) / determinant

        return false, Point2D(x, y)
    end
end

function order_intersection_points(track::Track, x1::Point2D, x2::Point2D)
    @unpack ϕ = track
    if isless(ϕ, π / 2)
        xi, xo = first(x1) < first(x2) ? (x1, x2) : (x2, x1)
    else
        xi, xo = first(x1) > first(x2) ? (x1, x2) : (x2, x1)
    end
    return xi, xo
end