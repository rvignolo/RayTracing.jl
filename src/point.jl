"""
    Point2D{T<:Real} <: FieldVector{2,T}

A 2D point in Euclidean space with coordinates (x, y).

This structure represents a point in two-dimensional space and inherits from
`FieldVector{2,T}`, providing vector arithmetic operations and StaticArrays
functionality for efficient geometric calculations in ray tracing.

## Fields
- `x::T`: x-coordinate
- `y::T`: y-coordinate

## Constructor
The constructor `Point2D(x, y)` automatically promotes the input arguments
to a common numeric type for type stability.

## Usage
```julia
# Create points with different numeric types
p1 = Point2D(1.0, 2.0)      # Float64
p2 = Point2D(1, 2)          # Promoted to common type
p3 = Point2D(1.0f0, 2.0f0)  # Float32

# Vector operations (inherited from FieldVector)
p4 = p1 + p2                # Addition
distance = norm(p1 - p2)    # Euclidean distance
```

Used extensively in ray tracing for track entry/exit points and segment endpoints.

See also: [`Track`](@ref), [`Segment`](@ref), [`advance_step`](@ref)
"""
struct Point2D{T<:Real} <: FieldVector{2,T}
    x::T
    y::T
end
Point2D(x, y) = Point2D(promote(x, y)...)

StaticArrays.similar_type(::Type{<:Point2D}, ::Type{T}, ::Size{(2,)}) where {T} = Point2D{T}
Base.convert(::Type{<:Point2D{T}}, arg::VectorValue{D}) where {D,T} = Point2D{T}(Tuple(arg))

@inline advance_step(x::Point2D, step::Real, ϕ::Real) = x + step * Point2D(cos(ϕ), sin(ϕ))