"""
    Vacuum
    Reflective
    Periodic

Boundary condition types for ray tracing simulations.

- `Vacuum`: Rays exit the domain without reflection
- `Reflective`: Rays are reflected back into the domain
- `Periodic`: Rays wrap around to the opposite boundary
"""
@enum BoundaryType begin
    Vacuum
    Reflective
    Periodic
end

"""
    BoundaryConditions{T<:BoundaryType}

Specifies boundary conditions for the four sides of a rectangular domain.

## Fields
- `top`: Boundary condition for the top edge
- `bottom`: Boundary condition for the bottom edge
- `right`: Boundary condition for the right edge
- `left`: Boundary condition for the left edge

## Examples
```julia
# All reflective boundaries
bcs = BoundaryConditions(top=Reflective, bottom=Reflective, right=Reflective, left=Reflective)

# Mixed boundary conditions
bcs = BoundaryConditions(top=Vacuum, bottom=Reflective, right=Periodic, left=Periodic)
```
"""
struct BoundaryConditions{T<:BoundaryType}
    top::T
    bottom::T
    right::T
    left::T
end

BoundaryConditions(; top=Vacuum, bottom=Vacuum, right=Vacuum, left=Vacuum) =
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
        error("Point do not lie in the boundary.")
    end

    return bc
end