
struct Mesh{M,K,B,N,C}
    model::M
    kdtree::K
    bounding_box::B
    node_cells::N
    cell_nodes::C
    # last_chosen_element # cache
end

function Mesh(model)
    grid = get_grid(model)
    kdtree = KDTree(grid)
    bounding_box = BoundingBox(grid)
    node_cells = get_faces(get_grid_topology(model), 0, num_cell_dims(model))
    cell_nodes = grid.cell_nodes
    return Mesh(model, kdtree, bounding_box, node_cells, cell_nodes)
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

struct BoundingBox{T}
    min::T
    max::T
end

"""
    BoundingBox(grid::UnstructuredGrid)

Finds the bounding box of a unstructured `grid`.
"""
function BoundingBox(grid::UnstructuredGrid{Dc,Dp,Tp}) where {Dc,Dp,Tp}
    nodes = get_node_coordinates(grid)

    # N = length(eltype(nodes)) # Dp
    # T = eltype(eltype(nodes)) # Tp

    xmin = MVector{Dp,Tp}(zeros(Dp)) # Dc or Dp?
    xmax = MVector{Dp,Tp}(zeros(Dp)) # Dc or Dp?

    for i in 1:Dp
        xs = getindex.(nodes, i)
        xmin[i] = min(xs...)
        xmax[i] = max(xs...)
    end

    bottomleft = convert(Point{Dp,Tp}, xmin)
    topright = convert(Point{Dp,Tp}, xmax)

    return BoundingBox(bottomleft, topright)
end

width(b::BoundingBox) = b.max[1] - b.min[1]
height(b::BoundingBox) = b.max[2] - b.min[2]

# retorna el elemento de la malla al que pertenece el punto en el espacio `x`
# mejor: function find_element(mesh, x) con mesh.kdtree y mesh.grid inside (also last_chosen_element)
# function find_element(grid, kdtree, )
function find_element(mesh, x, factor=2)
    @unpack model, kdtree, node_cells, cell_nodes = mesh
    grid = get_grid(model)

    # 1. test last chosen element (look cache)
    # element = mesh.last_chosen_element
    # if isassigned(element) end

    # ahora si arrancamos
    # 1. find nearest node
    nn_id, nn_distance = nn(kdtree, x)

    # 2. dada el id del nearest node, tomamos las cells que tienen ese nodo
    # TODO: se puede hacer esto sin allocar?
    cells = node_cells[nn_id]

    # mirando cada cell index, checkequemos si alguna contiene al punto en el espacio
    for cell in cells
        if point_in_element(grid, cell_nodes[cell], x)
            return cell
        end
    end

    # si la malla esta deformada, i.e. los elementos que pertenecen al nodo mas cercano no
    # contienen al punto, tenemos que ampliar la busqueda
    # TODO: en lugar de incrementar el radio de busqueda, podriamos usar knn, y recorrer
    # mas nodos. Esto me parece una mejor idea...
    nns_ids = inrange(kdtree, x, factor * nn_distance, true)
    for node_id in nns_ids
        cells = node_cells[node_id]
        for cell in cells
            if point_in_element(grid, cell_nodes[cell], x)
                return cell
            end
        end
    end

    return -1
end

# por ahora asumimos que son triangles siempre
point_in_element(grid, cell_nodes, x) = point_in_triangle(grid, cell_nodes, x)

function point_in_triangle(grid, nodes, x)

    node_coordinates = get_node_coordinates(grid)

    x1, y1 = node_coordinates[nodes[1]]
    x2, y2 = node_coordinates[nodes[2]]
    x3, y3 = node_coordinates[nodes[3]]

    # TODO: this is a first approach, use linear LinearAlgebra.jl or LazySets.jl later
    λ1 = ((y2 - y3) * (x[1] - x3) + (x3 - x2 ) * (x[2] - y3)) / ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1- y3))
    λ2 = ((y3 - y1) * (x[1] - x3) + (x1 - x3 ) * (x[2] - y3)) / ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1- y3))
    λ3 = 1 - λ1 - λ2

    # devuelve true si cae dentro o en la superficie del triangulo
    T = typeof(λ1)
    tol = sqrt(eps(T)) # por ahi es muy chica... ver de incrementarla si es necesario
    domain = ClosedInterval{T}(zero(T) - tol, one(T) + tol)
    # domain = OpenInterval{T}(zero(T) - tol, one(T) + tol) # devuelve true solo si es dentro
    return λ1 in domain && λ2 in domain && λ3 in domain
end