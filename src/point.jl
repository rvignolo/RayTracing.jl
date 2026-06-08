"""
    Point2D{T<:Real} <: FieldVector{2,T}

A 2D point in Euclidean space with coordinates (x, y).

This structure represents a point in two-dimensional space and inherits from
`FieldVector{2,T}`, providing vector arithmetic operations and StaticArrays functionality
for efficient geometric calculations in ray tracing.

## Type Parameters
- `T<:Real`: The numeric type for coordinates (Float64, Float32, etc.)

## Fields
- `x::T`: x-coordinate
- `y::T`: y-coordinate

## Constructors

### Basic Constructor
```julia
Point2D(x, y)
```
Automatically promotes input arguments to a common numeric type for type stability.

### Type-Specific Constructor
```julia
Point2D{T}(x, y)
```
Creates a point with specific numeric type T.

## Geometric Operations

The type inherits all vector operations from `FieldVector`:
- **Arithmetic**: `+`, `-`, `*`, `/`
- **Scalar operations**: `*`, `/`
- **Comparison**: `==`, `≈`, `<`, `≤`, `>`, `≥`
- **Norms**: `norm`, `norm_sqr`
- **Distance**: `distance`, `distance_sqr`

## Usage Examples

```julia
# Create points with different numeric types
p1 = Point2D(1.0, 2.0)      # Float64
p2 = Point2D(1, 2)          # Promoted to common type
p3 = Point2D{Float32}(1.0f0, 2.0f0)  # Float32

# Vector operations
p4 = p1 + p2                # Addition
p5 = p1 - p2                # Subtraction
p6 = 2.0 * p1               # Scalar multiplication
distance = norm(p1 - p2)    # Euclidean distance

# Geometric operations
midpoint = (p1 + p2) / 2    # Midpoint between two points
angle = atan(p1.y, p1.x)    # Angle from origin
```

## Performance Characteristics

- **Zero allocation**: All operations are allocation-free
- **SIMD friendly**: StaticArrays enable vectorization
- **Type stable**: Consistent return types across operations
- **Inlined**: Critical operations are marked with `@inline`

## Notes

- Used extensively in ray tracing for track entry/exit points and segment endpoints
- Compatible with Gridap's `VectorValue` type for seamless integration
- Supports automatic type promotion for mixed numeric types
- All geometric operations are numerically robust

## See Also

- [`Track`](@ref): Ray trajectory with entry/exit points
- [`Segment`](@ref): Line segment with endpoints
- [`advance_step`](@ref): Move point along direction
- [`distance`](@ref): Euclidean distance between points
- [`midpoint`](@ref): Midpoint between two points
- [`angle`](@ref): Angle between points or vectors
"""
struct Point2D{T<:Real} <: FieldVector{2,T}
    x::T
    y::T
end

# Constructors
Point2D(x, y) = Point2D(promote(x, y)...)

# Type conversion and compatibility
StaticArrays.similar_type(::Type{<:Point2D}, ::Type{T}, ::Size{(2,)}) where {T} = Point2D{T}
Base.convert(::Type{<:Point2D{T}}, arg::VectorValue{D}) where {D,T} = Point2D{T}(Tuple(arg))

# Gridap compatibility
Base.convert(::Type{VectorValue{2,T}}, p::Point2D{T}) where {T} = VectorValue{2,T}(p.x, p.y)

"""
    advance_step(x::Point2D, step::Real, ϕ::Real) -> Point2D

Advances a point by a given step size in the direction specified by angle ϕ.

## Arguments
- `x::Point2D`: Starting point
- `step::Real`: Distance to advance (positive for forward, negative for backward)
- `ϕ::Real`: Direction angle in radians (0 = positive x-axis, π/2 = positive y-axis)

## Returns
- `Point2D`: New point after advancing

## Algorithm
Computes the new position using vector addition:
```
x_new = x + step * [cos(ϕ), sin(ϕ)]
```

## Examples
```julia
p = Point2D(0.0, 0.0)
p_forward = advance_step(p, 1.0, 0.0)      # Point2D(1.0, 0.0)
p_up = advance_step(p, 1.0, π/2)           # Point2D(0.0, 1.0)
p_backward = advance_step(p, -1.0, 0.0)    # Point2D(-1.0, 0.0)
```

## Notes
- The angle ϕ is in radians (not degrees)
- Positive step moves in the direction of ϕ
- Negative step moves in the opposite direction
- Uses efficient trigonometric functions
"""
@inline advance_step(x::Point2D, step::Real, ϕ::Real) = x + step * Point2D(cos(ϕ), sin(ϕ))

"""
    distance(p1::Point2D, p2::Point2D) -> Real

Computes the Euclidean distance between two points.

## Arguments
- `p1::Point2D`: First point
- `p2::Point2D`: Second point

## Returns
- `Real`: Euclidean distance between the points

## Algorithm
Uses the Euclidean distance formula: √((x₂-x₁)² + (y₂-y₁)²)

## Examples
```julia
p1 = Point2D(0.0, 0.0)
p2 = Point2D(3.0, 4.0)
d = distance(p1, p2)  # 5.0 (3-4-5 triangle)
```

## Performance
- Uses `norm` from StaticArrays for optimal performance
- Avoids square root when possible via `distance_sqr`
"""
@inline distance(p1::Point2D, p2::Point2D) = norm(p1 - p2)

"""
    distance_sqr(p1::Point2D, p2::Point2D) -> Real

Computes the squared Euclidean distance between two points.

## Arguments
- `p1::Point2D`: First point
- `p2::Point2D`: Second point

## Returns
- `Real`: Squared Euclidean distance between the points

## Algorithm
Uses the squared distance formula: (x₂-x₁)² + (y₂-y₁)²

## Examples
```julia
p1 = Point2D(0.0, 0.0)
p2 = Point2D(3.0, 4.0)
d_sqr = distance_sqr(p1, p2)  # 25.0 (avoiding square root)
```

## Performance
- Faster than `distance` when exact distance is not needed
- Useful for distance comparisons (avoiding square root)
- Used in nearest neighbor searches and spatial queries
"""
@inline distance_sqr(p1::Point2D, p2::Point2D) = norm_sqr(p1 - p2)

"""
    midpoint(p1::Point2D, p2::Point2D) -> Point2D

Computes the midpoint between two points.

## Arguments
- `p1::Point2D`: First point
- `p2::Point2D`: Second point

## Returns
- `Point2D`: Midpoint between p1 and p2

## Algorithm
Computes the arithmetic mean: (p1 + p2) / 2

## Examples
```julia
p1 = Point2D(0.0, 0.0)
p2 = Point2D(4.0, 6.0)
mid = midpoint(p1, p2)  # Point2D(2.0, 3.0)
```

## Notes
- The midpoint is equidistant from both input points
- Used in geometric algorithms and mesh operations
"""
@inline midpoint(p1::Point2D, p2::Point2D) = (p1 + p2) / 2

"""
    angle(p1::Point2D, p2::Point2D) -> Real

Computes the angle between two points relative to the positive x-axis.

## Arguments
- `p1::Point2D`: First point (origin)
- `p2::Point2D`: Second point

## Returns
- `Real`: Angle in radians from p1 to p2

## Algorithm
Uses `atan` function: atan(y₂-y₁, x₂-x₁)

## Examples
```julia
p1 = Point2D(0.0, 0.0)
p2 = Point2D(1.0, 1.0)
ϕ = angle(p1, p2)  # π/4 radians (45 degrees)
```

## Notes
- Returns angle in radians (not degrees)
- Range: [-π, π]
- 0 = positive x-axis, π/2 = positive y-axis
- Used in ray tracing for azimuthal angle calculations
"""
@inline angle(p1::Point2D, p2::Point2D) = atan(p2.y - p1.y, p2.x - p1.x)

"""
    angle_from_origin(p::Point2D) -> Real

Computes the angle of a point relative to the origin.

## Arguments
- `p::Point2D`: Point to compute angle for

## Returns
- `Real`: Angle in radians from origin to point

## Examples
```julia
p = Point2D(1.0, 1.0)
ϕ = angle_from_origin(p)  # π/4 radians (45 degrees)
```

## Notes
- Equivalent to `angle(Point2D(0,0), p)`
- Used in polar coordinate conversions
"""
@inline angle_from_origin(p::Point2D) = atan(p.y, p.x)

"""
    rotate(p::Point2D, ϕ::Real) -> Point2D

Rotates a point around the origin by angle ϕ.

## Arguments
- `p::Point2D`: Point to rotate
- `ϕ::Real`: Rotation angle in radians (positive = counter-clockwise)

## Returns
- `Point2D`: Rotated point

## Algorithm
Uses rotation matrix: [cos(ϕ) -sin(ϕ); sin(ϕ) cos(ϕ)] * [x; y]

## Examples
```julia
p = Point2D(1.0, 0.0)
p_rotated = rotate(p, π/2)  # Point2D(0.0, 1.0)
```

## Notes
- Positive angle rotates counter-clockwise
- Used in geometric transformations and coordinate system rotations
"""
@inline function rotate(p::Point2D, ϕ::Real)
    c, s = cos(ϕ), sin(ϕ)
    return Point2D(c * p.x - s * p.y, s * p.x + c * p.y)
end

"""
    rotate_around(p::Point2D, center::Point2D, ϕ::Real) -> Point2D

Rotates a point around a specified center by angle ϕ.

## Arguments
- `p::Point2D`: Point to rotate
- `center::Point2D`: Center of rotation
- `ϕ::Real`: Rotation angle in radians

## Returns
- `Point2D`: Rotated point

## Algorithm
1. Translate point relative to center
2. Rotate around origin
3. Translate back

## Examples
```julia
p = Point2D(2.0, 0.0)
center = Point2D(1.0, 0.0)
p_rotated = rotate_around(p, center, π/2)  # Point2D(1.0, 1.0)
```

## Notes
- Used in geometric transformations and mesh operations
"""
@inline function rotate_around(p::Point2D, center::Point2D, ϕ::Real)
    translated = p - center
    rotated = rotate(translated, ϕ)
    return rotated + center
end

"""
    is_approx(p1::Point2D, p2::Point2D; atol::Real=0, rtol::Real=Base.rtoldefault(...)) -> Bool

Checks if two points are approximately equal within specified tolerances.

## Arguments
- `p1::Point2D`: First point
- `p2::Point2D`: Second point
- `atol::Real=0`: Absolute tolerance
- `rtol`: Relative tolerance. By default, uses the same coordinate-type default as
  `Base.isapprox`.

## Returns
- `Bool`: `true` if points are approximately equal

## Algorithm
Uses `isapprox` for each coordinate component

## Examples
```julia
p1 = Point2D(1.0, 2.0)
p2 = Point2D(1.0 + 1e-10, 2.0 + 1e-10)
is_close = is_approx(p1, p2)  # true
```

## Notes
- Useful for numerical comparisons in geometric algorithms
- Handles floating-point precision issues
"""
@inline function is_approx(
    p1::Point2D,
    p2::Point2D;
    atol::Real=0,
    rtol::Real=Base.rtoldefault(p1.x, p2.x, atol)
)
    return isapprox(p1.x, p2.x, atol=atol, rtol=rtol) &&
           isapprox(p1.y, p2.y, atol=atol, rtol=rtol)
end

# Additional utility functions for ray tracing

"""
    direction_vector(ϕ::Real) -> Point2D

Creates a unit vector in the direction specified by angle ϕ.

## Arguments
- `ϕ::Real`: Direction angle in radians

## Returns
- `Point2D`: Unit vector [cos(ϕ), sin(ϕ)]

## Examples
```julia
dir = direction_vector(0.0)      # Point2D(1.0, 0.0)
dir = direction_vector(π/2)      # Point2D(0.0, 1.0)
dir = direction_vector(π)        # Point2D(-1.0, 0.0)
```

## Notes
- Used in ray tracing for direction vectors
- Always has unit length (norm = 1)
"""
@inline direction_vector(ϕ::Real) = Point2D(cos(ϕ), sin(ϕ))

"""
    perpendicular(p::Point2D) -> Point2D

Returns a vector perpendicular to the given point/vector.

## Arguments
- `p::Point2D`: Input point/vector

## Returns
- `Point2D`: Perpendicular vector [-y, x]

## Examples
```julia
p = Point2D(1.0, 0.0)
perp = perpendicular(p)  # Point2D(0.0, 1.0)
```

## Notes
- Returns the 90-degree counter-clockwise rotation
- Used in geometric algorithms and normal vector calculations
"""
@inline perpendicular(p::Point2D) = Point2D(-p.y, p.x)
