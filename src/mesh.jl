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
- `node_cells::N`: Mapping from node indices to the set of cells containing each node
- `cell_nodes::C`: Mapping from cell indices to the ordered list of node indices defining
  each cell
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
# Create a mesh from a geometric model from Gridap
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
    width(mesh::Mesh)

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
    node_cells = get_faces(get_grid_topology(model), 0, num_cell_dims(model))
    cell_nodes = get_cell_node_ids(grid)
    bb_min, bb_max = bounding_box(grid)
    return Mesh(model, kdtree, node_cells, cell_nodes, bb_min, bb_max)
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

    xmin = MVector{Dp,Tp}(zeros(Tp, Dp))
    xmax = MVector{Dp,Tp}(zeros(Tp, Dp))

    for i in 1:Dp
        xs = getindex.(nodes, i)
        xmin[i] = min(xs...)
        xmax[i] = max(xs...)
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
    return any(isapprox.(p, bb_min, atol=atol)) || any(isapprox.(p, bb_max, atol=atol))
end

# Use dispatch once I get the info about the element type using Gridap topology.
point_in_element(mesh::Mesh, node_ids::AbstractVector{<:Int32}, x::Point2D) =
    point_in_triangle(mesh, node_ids, x)

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
# Find which element contains a specific point
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

    # get the nearest node id closest to `x`
    nn_id, _ = nn(kdtree, x)

    # get the associated cell ids that contain the nearest node
    cell_ids = node_cells[nn_id]

    # loop over those cells until the element containing `x` is found
    element = _find_element_in_cells(mesh, cell_ids, x)
    element != -one(eltype(cell_ids)) && return element

    # the mesh might be deformed, i.e. the cells that contain the nearest node do not
    # contain the point `x`. In that case, we need to search for more nodes.
    nn_ids, _ = knn(kdtree, x, k, true, i -> isequal(i, nn_id)) # TODO: cache!
    for node_id in nn_ids
        cell_ids = node_cells[node_id]
        element = _find_element_in_cells(mesh, cell_ids, x)
        element != -one(eltype(cell_ids)) && return element
    end

    return -one(eltype(cell_ids))
end

"""
    point_in_triangle(mesh::Mesh, node_ids, x) -> Bool

Checks whether a given point `x` lies inside, the edge or corner of the triangle given by
its node coordinates ids `node_ids`.
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

    # return true if it lies in or on the triangle
    T = eltype(λ)
    tol = sqrt(eps(T))
    domain = ClosedInterval{T}(zero(T) - tol, one(T) + tol)
    return λ[1] in domain && λ[2] in domain && λ[3] in domain
    # return all(in.(λ, Ref(domain))) # allocates
end

"""
    point_in_quadrangle(mesh::Mesh, node_ids, x) -> Bool

Checks whether a given point `x` lies inside the quadrangle given by its node coordinates
ids `node_ids`.
"""
function point_in_quadrangle(mesh::Mesh, node_ids::AbstractVector{<:Int32}, x::Point2D)

    # triangle node ids
    t_node_ids = @MVector zeros(3)

    # look on 4 triangles because we do not know the order of the nodes
    for i in 1:4
        for j in 1:3
            k = mod1(i + j - 1, 4) # k = (i + j - 2) % 4 + 1, k = mod(i + j - 1, 1:4)
            t_node_ids[j] = node_ids[k]
        end
        if point_in_triangle(mesh, t_node_ids, x)
            return true
        end
    end

    return false
end