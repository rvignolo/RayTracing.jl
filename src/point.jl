"""
    Point2D{T<:Real} <: FieldVector{2,T}

Representation of a point in two-dimensional Euclidean space.
"""
struct Point2D{T<:Real} <: FieldVector{2,T}
    x::T
    y::T
end
Point2D(x, y) = Point2D(promote(x, y)...)

StaticArrays.similar_type(::Type{<:Point2D}, ::Type{T}, ::Size{(2,)}) where {T} = Point2D{T}
Base.convert(::Type{<:Point2D{T}}, arg::VectorValue{D}) where {D,T} = Point2D{T}(Tuple(arg))

@inline advance_step(x::Point2D, step::Real, ϕ::Real) = x + step * Point2D(cos(ϕ), sin(ϕ))