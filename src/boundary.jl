"""
    Vacuum
    Reflective
    Periodic

Boundary condition types for ray tracing simulations in neutron transport.

These boundary conditions define how rays interact with the domain boundaries during ray tracing
calculations.

## Boundary Types

- **`Vacuum`**: Rays exit the domain without reflection (absorbing boundary)
- **`Reflective`**: Rays are reflected back into the domain (specular reflection)
- **`Periodic`**: Rays wrap around to the opposite boundary (periodic domain)

## Physical Interpretation

- **Vacuum**: Represents an absorbing boundary where neutrons escape the system
- **Reflective**: Represents a perfectly reflecting surface (e.g., symmetry planes)
- **Periodic**: Represents a repeating domain structure (e.g., lattice cells)
"""
@enum BoundaryType begin
    Vacuum
    Reflective
    Periodic
end

"""
    BoundaryConditions{T<:BoundaryType}

Specifies boundary conditions for the four sides of a rectangular domain.

This structure defines how rays interact with each boundary of the computational domain during ray
tracing calculations. Each side can have different boundary conditions to model various physical
scenarios.

## Type Parameters
- `T<:BoundaryType`: Type of boundary condition (Vacuum, Reflective, or Periodic)

## Fields
- `top::T`: Boundary condition for the top edge (y = y_max)
- `bottom::T`: Boundary condition for the bottom edge (y = y_min)
- `right::T`: Boundary condition for the right edge (x = x_max)
- `left::T`: Boundary condition for the left edge (x = x_min)

## Common Configurations

### All Reflective (Symmetric Domain)
```julia
bcs = BoundaryConditions(top=Reflective, bottom=Reflective, right=Reflective, left=Reflective)
```
Useful for symmetric problems where the domain represents a quarter or half of the full geometry.

### All Vacuum (Open Domain)
```julia
bcs = BoundaryConditions(top=Vacuum, bottom=Vacuum, right=Vacuum, left=Vacuum)
```
Represents an open domain where rays can escape in all directions.

### Mixed Conditions (Realistic Scenarios)
```julia
# Symmetry on left/right, open on top/bottom
bcs = BoundaryConditions(top=Vacuum, bottom=Vacuum, right=Reflective, left=Reflective)

# Periodic in x-direction, reflective in y-direction
bcs = BoundaryConditions(top=Reflective, bottom=Reflective, right=Periodic, left=Periodic)
```

## Notes

- **Reflective boundaries** are commonly used to exploit geometric symmetry
- **Periodic boundaries** are useful for lattice problems and infinite domains
- **Vacuum boundaries** model realistic scenarios where neutrons can escape
- All sides can have different conditions to model complex geometries

See also: [`boundary_condition`](@ref), [`BoundaryType`](@ref)
"""
struct BoundaryConditions{T<:BoundaryType}
    top::T
    bottom::T
    right::T
    left::T
end

"""
    BoundaryConditions(; top=Vacuum, bottom=Vacuum, right=Vacuum, left=Vacuum)

Construct boundary conditions with default vacuum boundaries on all sides.

## Arguments
- `top`: Boundary condition for top edge (default: Vacuum)
- `bottom`: Boundary condition for bottom edge (default: Vacuum)
- `right`: Boundary condition for right edge (default: Vacuum)
- `left`: Boundary condition for left edge (default: Vacuum)

## Returns
- `BoundaryConditions` structure with specified boundary conditions

## Examples

```julia
# Default (all vacuum)
bcs = BoundaryConditions()

# All reflective
bcs = BoundaryConditions(top=Reflective, bottom=Reflective, right=Reflective, left=Reflective)

# Mixed conditions
bcs = BoundaryConditions(top=Vacuum, bottom=Reflective, right=Periodic, left=Periodic)
```
"""
BoundaryConditions(; top=Vacuum, bottom=Vacuum, right=Vacuum, left=Vacuum) =
    BoundaryConditions(top, bottom, right, left)

"""
    boundary_condition(x::Point2D, sides, bcs::BoundaryConditions)

Determine the boundary condition at a given point on the domain boundary.

This function identifies which boundary a point lies on and returns the corresponding boundary
condition. It's used during ray tracing to determine how rays interact with domain boundaries.

## Arguments
- `x::Point2D`: Point on the domain boundary
- `sides`: Named tuple containing boundary segments (top, bottom, right, left)
- `bcs::BoundaryConditions`: Boundary conditions for all sides

## Returns
- `BoundaryType`: The boundary condition at the given point

## Throws
- `ArgumentError`: If the point does not lie on any boundary

## Example

```julia
# Define boundary conditions
bcs = BoundaryConditions(top=Vacuum, bottom=Reflective, right=Periodic, left=Periodic)

# Define domain sides
sides = (top=Segment(p1, p2), bottom=Segment(p3, p4),
         right=Segment(p2, p5), left=Segment(p4, p1))

# Check boundary condition at a point
bc = boundary_condition(point, sides, bcs)
```

## Notes

- The point must lie exactly on one of the boundary segments
- The function checks boundaries in order: top, bottom, right, left
- Returns the first matching boundary condition
"""
function boundary_condition(x::Point2D, sides, bcs::BoundaryConditions)
    if x in sides.top
        return bcs.top
    elseif x in sides.bottom
        return bcs.bottom
    elseif x in sides.right
        return bcs.right
    elseif x in sides.left
        return bcs.left
    else
        throw(ArgumentError("Point $(x) does not lie on any boundary. " *
                            "Ensure the point is exactly on one of the domain edges."))
    end
end

# Convenience constructors for common boundary configurations
"""
    reflective_boundaries()

Create boundary conditions with all sides set to reflective.

Useful for symmetric problems where the domain represents a portion of a larger symmetric geometry.

## Returns
- `BoundaryConditions` with all sides set to `Reflective`

## Example
```julia
bcs = reflective_boundaries()
# Equivalent to: BoundaryConditions(top=Reflective, bottom=Reflective, right=Reflective, left=Reflective)
```
"""
reflective_boundaries() = BoundaryConditions(top=Reflective, bottom=Reflective, right=Reflective, left=Reflective)

"""
    vacuum_boundaries()

Create boundary conditions with all sides set to vacuum.

Represents an open domain where rays can escape in all directions.

## Returns
- `BoundaryConditions` with all sides set to `Vacuum`

## Example
```julia
bcs = vacuum_boundaries()
# Equivalent to: BoundaryConditions() (default constructor)
```
"""
vacuum_boundaries() = BoundaryConditions()

"""
    periodic_boundaries()

Create boundary conditions with all sides set to periodic.

Useful for infinite domain problems or lattice structures.

## Returns
- `BoundaryConditions` with all sides set to `Periodic`

## Example
```julia
bcs = periodic_boundaries()
# Equivalent to: BoundaryConditions(top=Periodic, bottom=Periodic, right=Periodic, left=Periodic)
```
"""
periodic_boundaries() = BoundaryConditions(top=Periodic, bottom=Periodic, right=Periodic, left=Periodic)