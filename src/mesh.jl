# TODO: use elements for any kind of object and use cell for the other!

"""
    Mesh

Holds mesh related attributes that are useful for the ray tracing algorithm, such as a
`kdtree` for nodes searches, the mesh bounding box and a mapping between nodes and cells
(`node_cells` and `cell_nodes`).
"""
struct Mesh{M,K,N,C,B}
    model::M
    kdtree::K
    node_cells::N
    cell_nodes::C
    bb_min::B
    bb_max::B
end

"""
    Mesh(model::UnstructuredDiscreteModel)

Builds a Mesh using data from a [`UnstructuredDiscreteModel`](@ref).
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
    KDTree(grid::UnstructuredGrid)

Builds a [`KDTree`](@ref) using node data from a `grid`.
"""
function KDTree(grid::UnstructuredGrid{Dc,Dp,Tp}) where {Dc,Dp,Tp}
    nodes = get_node_coordinates(grid)
    snodes = convert.(SVector{Dp,Tp}, nodes) # Dc or Dp?
    return KDTree(snodes)
end

num_dims(mesh::Mesh) = num_dims(mesh.model)
num_cells(mesh::Mesh) = num_cells(mesh.model)
num_nodes(mesh::Mesh) = num_nodes(mesh.model)

"""
    bounding_box(grid::UnstructuredGrid)

Returns the bounding box of an unstructured `grid`.
"""
function bounding_box(grid::UnstructuredGrid{Dc,Dp,Tp}) where {Dc,Dp,Tp}
    nodes = get_node_coordinates(grid)

    xmin = MVector{Dp,Tp}(zeros(Tp, Dp)) # Dc or Dp?
    xmax = MVector{Dp,Tp}(zeros(Tp, Dp))

    for i in 1:Dp
        xs = getindex.(nodes, i)
        xmin[i] = min(xs...)
        xmax[i] = max(xs...)
    end

    bottomleft = convert(Point2D{Tp}, xmin)
    topright = convert(Point2D{Tp}, xmax)

    return bottomleft, topright
end

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
    inboundary(mesh::Mesh, x::Point2D, [atol::Real=0]) -> Bool

Checks if a point `x` lies in the boundary of the `mesh` with certain absolute tolerance
`atol`.
"""
function inboundary(mesh::Mesh, x::Point2D, atol::Real=0)
    @unpack bb_min, bb_max = mesh
    return isapprox(x[1], bb_max[1], atol=atol) || isapprox(x[1], bb_min[1], atol=atol) ||
           isapprox(x[2], bb_max[2], atol=atol) || isapprox(x[2], bb_min[2], atol=atol)
end

"""
    find_element(mesh::Mesh, x::Point, k::Int=2)

Finds the cell or element of the `mesh` that contains a given point `x` by searching nodes
and elements using nearest neighbor searches.
"""
function find_element(mesh::Mesh, x::Point2D, k::Int=2)
    @unpack model, kdtree, node_cells, cell_nodes = mesh

    # get the nearest node id closest to `x`
    nn_id, nn_distance = nn(kdtree, x)

    # get the associated cell ids that contain the nearest node
    cell_ids = node_cells[nn_id]

    # loop over those cells until we find the one that contains the point `x`
    for cell_id in cell_ids
        node_ids = cell_nodes[cell_id]
        if point_in_element(mesh, node_ids, x)
            return cell_id
        end
    end

    # the mesh might be deformed, i.e. the cells that contain the nearest node do not
    # contain the point `x`. In that case, we need to search for more nodes. We could either
    # use `inrange` or `knn` for this purpose.
    nn_ids, nn_distances = knn(kdtree, x, k, true, i -> isequal(i, nn_id)) # TODO: cache!
    for node_id in nn_ids
        cell_ids = node_cells[node_id]
        for cell_id in cell_ids
            node_ids = cell_nodes[cell_id]
            if point_in_element(mesh, node_ids, x)
                return cell_id
            end
        end
    end

    # nn_ids = inrange(kdtree, x, factor * nn_distance, true)
    # for node_id in nn_ids
    #     cell_ids = node_cells[node_id]
    #     for cell_id in cell_ids
    #         node_ids = cell_nodes[cell_id]
    #         if point_in_element(grid, node_ids, x)
    #             return cell_id
    #         end
    #     end
    # end

    return -1
end

# Use dispatch once I get the info about the element type using Gridap topology.
point_in_element(mesh::Mesh, node_ids::AbstractVector{<:Int32}, x::Point2D) =
    point_in_triangle(mesh, node_ids, x)

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