struct BoundingBox{I,A}
    min::I
    max::A
end

# function BoundingBox(model::GridapGmsh.UnstructuredDiscreteModel)
function BoundingBox(model)
    grid = get_grid(model)
    nodes = grid.node_coordinates

    xmin = ymin = +Inf
    xmax = ymax = -Inf
    for node in nodes
        xmin = node[1] < xmin ? node[1] : xmin
        ymin = node[2] < ymin ? node[2] : ymin
        xmax = node[1] > xmax ? node[1] : xmax
        ymax = node[2] > ymax ? node[2] : ymax
    end

    T = eltype(eltype(nodes))

    return BoundingBox(Point{2,T}(xmin, ymin), Point{2,T}(xmax, ymax))
end

width(b::BoundingBox) = b.max[1] - b.min[1]
height(b::BoundingBox) = b.max[2] - b.min[2]