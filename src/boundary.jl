@enum BoundaryType begin
    Vaccum
    Reflective
    Periodic
end

struct BoundaryConditions{T<:BoundaryType}
    top::T
    bottom::T
    right::T
    left::T
end

BoundaryConditions(; top=Vaccum, bottom=Vaccum, right=Vaccum, left=Vaccum) =
    BoundaryConditions(top, bottom, right, left)

function boundary_condition(x::Point2D, sides, bcs::BoundaryConditions)

    if x in sides.top
        bc = bcs.top
    elseif x in sides.bottom
        bc = bcs.bottom
    elseif x in sides.right
        bc = bcs.right
    elseif x in sides.left
        bc = bcs.left
    else
        error("Entry or exit point do not lie in the mesh boundary.")
    end

    return bc
end