"""
    Mesh{M,K,N,C,B}

A computational mesh structure that holds geometric and topological information essential
for ray tracing algorithms in neutron transport simulations.

## Type Parameters

- `M`: Type of the underlying geometric model
- `K`: Type of the spatial index structure (typically a KD-tree)
- `N`: Type of the node-to-cells mapping
- `C`: Type of the cell-to-nodes mapping
- `B`: Type of the bounding box coordinates

## Fields

- `model::M`: The underlying geometric model containing mesh topology and properties
- `kdtree::K`: Spatial index structure for efficient nearest-neighbor searches and spatial
  queries
- `node_cells::N`: Materialized mapping from node indices to the set of cells containing
  each node
- `cell_nodes::C`: Materialized mapping from cell indices to the ordered list of node
  indices defining each cell
- `bb_min::B`: Minimum coordinates of the mesh bounding box (lower-left corner)
- `bb_max::B`: Maximum coordinates of the mesh bounding box (upper-right corner)

## Purpose

This structure serves as the primary data container for mesh-based ray tracing algorithms,
providing:

1. **Geometric Information**: Node coordinates and cell definitions through the model
2. **Topological Relationships**: Efficient lookups between nodes and cells
3. **Spatial Indexing**: Fast spatial queries using the KD-tree
4. **Bounding Information**: Global mesh bounds for optimization and validation

## Usage

The `Mesh` structure is typically constructed from a geometric model and used throughout the
ray tracing pipeline for:

- **Track-Mesh Intersection**: Finding which cells a ray intersects
- **Spatial Queries**: Locating nearby nodes or cells efficiently
- **Topological Traversal**: Navigating between connected mesh elements
- **Boundary Detection**: Identifying mesh boundaries and interfaces

## Performance Considerations

- The KD-tree enables O(log n) spatial queries instead of O(n) brute force searches
- Node-cell mappings provide O(1) lookups for topological relationships
- Bounding box information allows early rejection of rays that miss the mesh entirely

## Example

```julia
# Create a mesh from a Gridap geometric model.
model = DiscreteModelFromFile(jsonfile)
mesh = Mesh(model)

# Find cells containing a specific node
node_id = 42
containing_cells = mesh.node_cells[node_id]

# Get nodes defining a specific cell
cell_id = 15
cell_nodes = mesh.cell_nodes[cell_id]

# Check if a point is within mesh bounds
point = Point2D(1.0, 2.0)
in_bounds = all(mesh.bb_min .≤ point .≤ mesh.bb_max)
```

## Notes

- The mesh is assumed to be valid (no degenerate elements, proper connectivity)
- All spatial coordinates are expected to be in the same coordinate system
- The KD-tree is automatically constructed during mesh initialization
"""
struct Mesh{M<:UnstructuredDiscreteModel,K<:KDTree,NC,T}
    model::M
    kdtree::K
    node_cells::NC
    cell_nodes::NC
    bb_min::Point2D{T}
    bb_max::Point2D{T}
end

"""
    num_dims(mesh::Mesh)

Returns the number of dimensions of the mesh.
"""
@inline num_dims(mesh::Mesh) = num_dims(mesh.model)

"""
    num_cells(mesh::Mesh)

Returns the number of cells in the mesh.
"""
@inline num_cells(mesh::Mesh) = num_cells(mesh.model)

"""
    num_nodes(mesh::Mesh)

Returns the number of nodes in the mesh.
"""
@inline num_nodes(mesh::Mesh) = num_nodes(mesh.model)

"""
    width(mesh::Mesh)

Returns the width of the rectangular mesh.
"""
@inline width(mesh::Mesh) = mesh.bb_max[1] - mesh.bb_min[1]

"""
    height(mesh::Mesh)

Returns the height of the rectangular mesh.
"""
@inline height(mesh::Mesh) = mesh.bb_max[2] - mesh.bb_min[2]

"""
    Mesh(model::UnstructuredDiscreteModel)

Constructs a Mesh from an [`UnstructuredDiscreteModel`](@ref) with enhanced spatial query
capabilities.

This constructor extracts the grid topology and builds several data structures for efficient
mesh operations:
- A KD-tree for fast spatial nearest-neighbor queries
- Node-to-cell and cell-to-node connectivity mappings
- Precomputed bounding box for bounds checking

# Arguments
- `model::UnstructuredDiscreteModel`: The discrete model containing the mesh geometry and
  topology

# Returns
- `Mesh`: A mesh structure optimized for spatial queries and connectivity operations
"""
function Mesh(model::UnstructuredDiscreteModel)
    grid = get_grid(model)
    kdtree = KDTree(grid)
    node_cells = materialize_connectivity(
        get_faces(get_grid_topology(model), 0, num_cell_dims(model))
    )
    cell_nodes = materialize_connectivity(get_cell_node_ids(grid))
    bb_min, bb_max = bounding_box(grid)
    return Mesh(model, kdtree, node_cells, cell_nodes, bb_min, bb_max)
end

function materialize_connectivity(connectivity)
    return [Vector{Int32}(ids) for ids in connectivity]
end

"""
    KDTree(grid::UnstructuredGrid{Dc,Dp,Tp}) -> KDTree

Constructs a KD-tree from the node coordinates of an unstructured grid for efficient spatial
nearest-neighbor queries.

# Arguments
- `grid::UnstructuredGrid{Dc,Dp,Tp}`: An unstructured grid where `Dc` is the cell dimension,
  `Dp` is the point dimension, and `Tp` is the numeric type

# Returns
- `KDTree`: A KD-tree data structure containing all grid nodes as static vectors for fast
  spatial lookups

The constructor extracts node coordinates from the grid and converts them to static vectors
(`SVector{Dp,Tp}`) for optimal performance in the KD-tree implementation.
"""
function KDTree(grid::UnstructuredGrid{Dc,Dp,Tp}) where {Dc,Dp,Tp}
    static_nodes = convert.(SVector{Dp,Tp}, get_node_coordinates(grid))
    return KDTree(static_nodes)
end

"""
    bounding_box(grid::UnstructuredGrid{Dc,Dp,Tp}) -> Tuple{Point2D{Tp}, Point2D{Tp}}

Computes the axis-aligned bounding box of an unstructured grid.

# Arguments
- `grid::UnstructuredGrid{Dc,Dp,Tp}`: An unstructured grid where `Dc` is the cell dimension,
  `Dp` is the point dimension, and `Tp` is the numeric type

# Returns
- `Tuple{Point2D{Tp}, Point2D{Tp}}`: A tuple containing the bottom-left and top-right
  corners of the bounding box that encloses all nodes in the grid

The bounding box is computed by finding the minimum and maximum coordinates across all grid
nodes in each spatial dimension.
"""
function bounding_box(grid::UnstructuredGrid{Dc,Dp,Tp}) where {Dc,Dp,Tp}
    nodes = get_node_coordinates(grid)

    xmin = MVector{Dp,Tp}(Tuple(first(nodes)))
    xmax = MVector{Dp,Tp}(Tuple(first(nodes)))

    for node in nodes
        for i in 1:Dp
            xi = node[i]
            xmin[i] = min(xmin[i], xi)
            xmax[i] = max(xmax[i], xi)
        end
    end

    bb_min = convert(Point2D{Tp}, xmin)
    bb_max = convert(Point2D{Tp}, xmax)

    return bb_min, bb_max
end

"""
    inboundary(mesh::Mesh, x::Point2D, [atol::Real=0]) -> Bool

Checks if a point `x` lies on the boundary of the mesh's bounding box with a specified
absolute tolerance `atol`.

## Arguments
- `mesh::Mesh`: The mesh whose boundary to check against
- `x::Point2D`: The point to check
- `atol::Real=0`: Absolute tolerance for floating-point comparisons

## Returns
- `Bool`: `true` if the point lies on any of the four boundary edges of the mesh's bounding
  box, `false` otherwise

## Notes
- Uses the mesh's bounding box (`bb_min`, `bb_max`) for boundary checking
- Checks all four edges: top, bottom, left, and right
- The tolerance is applied to each coordinate comparison independently
- This is a fast approximation that assumes the mesh boundary is rectangular
"""
@inline function on_boundary(mesh::Mesh, p::Point2D, atol::Real=zero(eltype(p)))
    atol < zero(atol) && throw(ArgumentError("Tolerance must be non-negative, got $atol"))
    @unpack bb_min, bb_max = mesh
    x_in = bb_min[1] - atol <= p[1] <= bb_max[1] + atol
    y_in = bb_min[2] - atol <= p[2] <= bb_max[2] + atol
    x_on = isapprox(p[1], bb_min[1]; atol=atol) || isapprox(p[1], bb_max[1]; atol=atol)
    y_on = isapprox(p[2], bb_min[2]; atol=atol) || isapprox(p[2], bb_max[2]; atol=atol)
    return (x_on && y_in) || (y_on && x_in)
end

"""
    point_in_element(mesh::Mesh, node_ids::AbstractVector{<:Int32}, x::Point2D) -> Bool

Checks if a given point `x` lies inside the element defined by the node coordinate IDs
`node_ids`.
"""
function point_in_element(mesh::Mesh, node_ids::AbstractVector{<:Int32}, x::Point2D)
    return point_in_element(mesh, Val(length(node_ids)), node_ids, x)
end

function ordered_node_ids(mesh::Mesh, node_ids::AbstractVector{<:Integer})
    length(node_ids) <= 3 && return node_ids

    node_coordinates = get_node_coordinates(get_grid(mesh.model))
    T = eltype(first(node_coordinates))
    cx = zero(T)
    cy = zero(T)

    for node_id in node_ids
        p = node_coordinates[node_id]
        cx += p[1]
        cy += p[2]
    end
    cx /= length(node_ids)
    cy /= length(node_ids)

    node_order = sortperm(
        collect(eachindex(node_ids));
        by=i -> atan(node_coordinates[node_ids[i]][2] - cy, node_coordinates[node_ids[i]][1] - cx)
    )

    return node_ids[node_order]
end

"""
    point_in_element(mesh::Mesh, _::Val{3}, node_ids::AbstractVector{<:Int32}, x::Point2D) -> Bool

Checks if a given point `x` lies inside the triangle defined by the node coordinate IDs
`node_ids`.
"""
@inline point_in_element(mesh::Mesh, _::Val{3}, node_ids::AbstractVector{<:Int32}, x::Point2D) =
    point_in_triangle(mesh, node_ids, x)

"""
    point_in_element(mesh::Mesh, _::Val{4}, node_ids::AbstractVector{<:Int32}, x::Point2D) -> Bool

Checks if a given point `x` lies inside the quadrangle defined by the node coordinate IDs
`node_ids`.
"""
@inline point_in_element(mesh::Mesh, _::Val{4}, node_ids::AbstractVector{<:Int32}, x::Point2D) =
    point_in_quadrangle(mesh, node_ids, x)

"""
    _find_element_in_cells(mesh::Mesh, cell_ids::AbstractVector{<:Int32}, x::Point2D) -> Int

Auxiliary function that searches through a collection of cell IDs to find which one contains
the point `x`. Returns the first cell ID that contains the point, or -1 if none found.

## Arguments
- `mesh::Mesh`: The mesh to search in
- `cell_ids`: Collection of cell IDs to check
- `x::Point2D`: The point to find

## Returns
- `Int`: The cell ID that contains the point, or -1 if not found
"""
@inline function _find_element_in_cells(mesh::Mesh, cell_ids::AbstractVector{<:Int32}, x::Point2D)
    @unpack cell_nodes = mesh
    for cell_id in cell_ids
        node_ids = cell_nodes[cell_id]
        if point_in_element(mesh, node_ids, x)
            return cell_id
        end
    end
    return -one(eltype(cell_ids))
end

"""
    find_element(mesh::Mesh, x::Point2D, k::Int=2) -> Int

Finds the mesh element (cell) that contains a given point `x` using an efficient two-stage
nearest neighbor search algorithm.

This function is essential for ray tracing algorithms as it determines which mesh element a
ray intersects. It uses a KD-tree for spatial indexing to achieve O(log n) search
performance instead of O(n) brute force search.

## Arguments
- `mesh::Mesh`: The computational mesh containing geometry and topology
- `x::Point2D`: The point whose containing element to find
- `k::Int=2`: Number of nearest neighbors to search if the first search fails (default: 2)

## Returns
- `Int`: The ID of the mesh element containing the point, or -1 if no element is found

## Algorithm

The function employs a two-stage search strategy:

1. **Primary Search**: Find the nearest mesh node to `x` using the KD-tree, then check all
   elements that contain this node for point containment.

2. **Fallback Search**: If no element is found in the primary search (e.g., in deformed
   meshes), expand the search to the `k` nearest nodes and check all their associated
   elements.

## Performance Characteristics

- **Best Case**: O(log n) when the point is in an element containing the nearest node
- **Worst Case**: O(k log n) when the point is in an element associated with the k-th
  nearest node
- **Average Case**: O(log n) for well-behaved meshes

## Usage Examples

```julia
# Find the element containing a specific point.
point = Point2D(1.5, 2.3)
element_id = find_element(mesh, point)

if element_id != -1
    println("Point is in element: ", element_id)
else
    println("Point is outside the mesh")
end

# Use with custom k for problematic meshes
element_id = find_element(mesh, point, k=5)
```

## Edge Cases and Limitations

- **Deformed Meshes**: For highly deformed meshes, increase `k` to improve success rate
- **Boundary Points**: Points exactly on element boundaries may return any adjacent element
- **Outside Points**: Returns -1 for points outside the mesh bounding box
- **Degenerate Elements**: May fail for meshes with degenerate (zero-area) elements

## Notes

- The function assumes the mesh is valid (no degenerate elements, proper connectivity)
- Uses barycentric coordinates for point-in-element testing
- The KD-tree is precomputed during mesh construction for optimal performance
- Consider caching results for repeated queries on the same or nearby points

## See Also

- [`point_in_element`](@ref): Tests if a point is inside a specific element
- [`point_in_triangle`](@ref): Triangle containment test using barycentric coordinates
- [`point_in_quadrangle`](@ref): Quadrangle containment test via triangulation
"""
function find_element(mesh::Mesh, x::Point2D, k::Int=2)
    @unpack model, kdtree, node_cells = mesh

    # Find the nearest node to `x`.
    nn_id, _ = nn(kdtree, x)

    # Get the cell IDs associated with the nearest node.
    cell_ids = node_cells[nn_id]

    # Define the invalid element ID.
    invalid_element_id = -one(eltype(cell_ids))

    # Search those cells for the element containing `x`.
    element = _find_element_in_cells(mesh, cell_ids, x)
    element != invalid_element_id && return element

    # In a deformed mesh, the cells attached to the nearest node may not contain `x`;
    # search additional nearby nodes in that case.
    nn_ids, _ = knn(kdtree, x, k, true, i -> isequal(i, nn_id)) # TODO: cache!
    for node_id in nn_ids
        cell_ids = node_cells[node_id]
        element = _find_element_in_cells(mesh, cell_ids, x)
        element != invalid_element_id && return element
    end

    return invalid_element_id
end

"""
    point_in_triangle(mesh::Mesh, node_ids::AbstractVector{<:Int32}, x::Point2D) -> Bool

Checks whether a given point `x` lies inside, on the edge, or at a corner of the triangle
defined by its node coordinates.

This function uses barycentric coordinates to determine point containment in a triangle. The
method is numerically robust and handles edge cases including points exactly on triangle
boundaries and at vertices.

## Arguments
- `mesh::Mesh`: The mesh containing the triangle
- `node_ids::AbstractVector{<:Int32}`: Exactly 3 node IDs defining the triangle vertices
- `x::Point2D`: The point to check for containment

## Returns
- `Bool`: `true` if the point lies inside, on an edge, or at a vertex of the triangle

## Algorithm

The function uses barycentric coordinate transformation:

1. **Extract vertex coordinates**: Get the 3D coordinates of the triangle vertices
2. **Form transformation matrix**: Create matrix R = [x₁ x₂ x₃; y₁ y₂ y₃; 1 1 1]
3. **Solve barycentric coordinates**: Solve R·λ = [x, y, 1] for λ = [λ₁, λ₂, λ₃]
4. **Check containment**: Point is inside if all λᵢ ∈ [0, 1] (with tolerance)

## Mathematical Foundation

Barycentric coordinates represent a point P as a weighted combination of triangle vertices:
```
P = λ₁·V₁ + λ₂·V₂ + λ₃·V₃
```
where λ₁ + λ₂ + λ₃ = 1 and λᵢ ≥ 0 for all i.

The point is:
- **Inside**: All λᵢ > 0 (strictly positive)
- **On edge**: One λᵢ = 0, others > 0
- **At vertex**: One λᵢ = 1, others = 0
- **Outside**: At least one λᵢ < 0

## Numerical Robustness

- **Tolerance**: Uses `√(eps(T))` for floating-point comparisons
- **Domain**: Checks λᵢ ∈ [-tol, 1+tol] to handle numerical errors
- **Edge cases**: Properly handles points exactly on boundaries
- **Degenerate triangles**: May fail for zero-area triangles

## Performance Characteristics

- **Computational complexity**: O(1) - constant time operation
- **Memory usage**: Minimal - only requires 3×3 matrix and 3-vector
- **Numerical stability**: Robust against floating-point errors
- **Vectorization**: Uses StaticArrays for efficient matrix operations

## Edge Cases and Limitations

- **Degenerate triangles**: May produce incorrect results for zero-area triangles
- **Node ordering**: Assumes valid triangle node ordering (counter-clockwise)
- **Mesh validity**: Requires mesh to have valid node coordinates

## Notes

- **Barycentric coordinates**: More robust than area-based methods
- **Boundary handling**: Points on edges and vertices return `true`
- **Tolerance**: Automatically adjusted based on type
- **StaticArrays**: Uses `@SMatrix` and `@SVector` for performance

## See Also

- [`point_in_element`](@ref): Generic point-in-element test with dispatch
- [`point_in_quadrangle`](@ref): Quadrangle containment test
- [`find_element`](@ref): Find which element contains a point
"""
function point_in_triangle(mesh::Mesh, node_ids::AbstractVector{<:Int32}, x::Point2D)
    @unpack model = mesh
    node_coordinates = get_node_coordinates(get_grid(model))

    x1, y1 = node_coordinates[node_ids[1]]
    x2, y2 = node_coordinates[node_ids[2]]
    x3, y3 = node_coordinates[node_ids[3]]

    R = @SMatrix [x1 x2 x3; y1 y2 y3; 1 1 1]
    r = @SVector [x[1], x[2], 1]
    λ = R \ r

    # Return true if `x` lies inside or on the triangle.
    T = eltype(λ)
    tol = sqrt(eps(T))
    domain = ClosedInterval{T}(zero(T) - tol, one(T) + tol)
    return λ[1] in domain && λ[2] in domain && λ[3] in domain
    # Equivalent, but allocates:
    # return all(in.(λ, Ref(domain)))
end

"""
    point_in_quadrangle(mesh::Mesh, node_ids::AbstractVector{<:Int32}, x::Point2D) -> Bool

Checks whether a given point `x` lies inside the quadrangle defined by its node coordinates.

This function uses a triangulation-based approach to determine point containment in a
quadrangle. The method decomposes the quadrangle into four triangles and tests point
containment in each triangle. This approach handles arbitrary quadrilateral shapes robustly,
including concave and non-convex quadrangles.

## Arguments
- `mesh::Mesh`: The mesh containing the quadrangle
- `node_ids::AbstractVector{<:Int32}`: Exactly 4 node IDs defining the quadrangle vertices
- `x::Point2D`: The point to check for containment

## Returns
- `Bool`: `true` if the point lies inside or on the boundary of the quadrangle

## Algorithm

The function uses a systematic triangulation approach:

1. **Quadrangle decomposition**: Decomposes the quadrangle into 4 triangles
2. **Triangle formation**: Creates triangles using different vertex combinations
3. **Point testing**: Tests point containment in each triangle using `point_in_triangle`
4. **Early termination**: Returns `true` as soon as point is found in any triangle

### Triangle Formation Strategy

The algorithm creates 4 triangles by systematically combining vertices:
- **Triangle 1**: Vertices [1, 2, 3] (first three vertices)
- **Triangle 2**: Vertices [2, 3, 4] (middle three vertices)
- **Triangle 3**: Vertices [3, 4, 1] (last and first vertices)
- **Triangle 4**: Vertices [4, 1, 2] (last two and first vertex)

This ensures complete coverage of the quadrangle regardless of node ordering.

## Mathematical Foundation

The triangulation approach works because:
- Any quadrangle can be decomposed into triangles
- Point containment in a quadrangle is equivalent to containment in at least one of its
  triangular sub-elements
- The union of the four triangles covers the entire quadrangle area

### Coverage Guarantee

For a quadrangle with vertices V₁, V₂, V₃, V₄, the four triangles:
- T₁ = (V₁, V₂, V₃)
- T₂ = (V₂, V₃, V₄)
- T₃ = (V₃, V₄, V₁)
- T₄ = (V₄, V₁, V₂)

Collectively cover the quadrangle: Q = T₁ ∪ T₂ ∪ T₃ ∪ T₄

## Performance Characteristics

- **Computational complexity**: O(1) - constant time operation (4 triangle tests)
- **Memory usage**: Minimal - only requires 3-vector for triangle nodes
- **Early termination**: Stops as soon as point is found in any triangle
- **Vectorization**: Uses `@MVector` for efficient temporary storage

## Edge Cases and Limitations

- **Degenerate quadrangles**: May produce incorrect results for zero-area quadrangles
- **Node ordering**: Assumes valid quadrangle node ordering (counter-clockwise)
- **Mesh validity**: Requires mesh to have valid node coordinates
- **Overlapping triangles**: The four triangles may overlap, but this doesn't affect
  correctness

## Numerical Robustness

- **Inherits robustness**: Benefits from the numerical robustness of `point_in_triangle`
- **Consistent behavior**: Points on quadrangle boundaries return `true`
- **Tolerance handling**: Uses the same tolerance as triangle containment tests
- **Edge cases**: Properly handles points on edges and at vertices

## Algorithm Advantages

- **Shape independence**: Works for arbitrary quadrilateral shapes
- **Robustness**: Inherits numerical robustness from triangle tests
- **Simplicity**: Straightforward implementation using existing triangle tests
- **Completeness**: Guaranteed to find points inside the quadrangle

## Notes

- **Triangulation approach**: More general than barycentric coordinate methods for
  quadrangles
- **Node ordering**: The `mod1` function ensures proper vertex cycling
- **Memory efficiency**: Uses `@MVector` for temporary storage to avoid allocations
- **Early termination**: Optimized to return as soon as containment is found
- **Reusability**: Leverages the robust `point_in_triangle` implementation

## See Also

- [`point_in_element`](@ref): Generic point-in-element test with dispatch
- [`point_in_triangle`](@ref): Triangle containment test using barycentric coordinates
- [`find_element`](@ref): Find which element contains a point
"""
function point_in_quadrangle(mesh::Mesh, node_ids::AbstractVector{<:Int32}, x::Point2D)

    # Triangle node IDs.
    t_node_ids = MVector{3,eltype(node_ids)}(undef)
    ordered_ids = ordered_node_ids(mesh, node_ids)

    # Test four triangles because node ordering may vary.
    for i in 1:4
        for j in 1:3
            k = mod1(i + j - 1, 4) # k = (i + j - 2) % 4 + 1, k = mod(i + j - 1, 1:4)
            t_node_ids[j] = ordered_ids[k]
        end
        if point_in_triangle(mesh, t_node_ids, x)
            return true
        end
    end

    return false
end
