"""
    Segment{T<:Real}

The intersection of a mesh element or cell and a [`Track`](@ref) yields to a `Segment` with
entry point `p`, exit point `q` and length `ℓ`.
"""
struct Segment{T<:Real}
    p::Point2D{T}
    q::Point2D{T}
    ℓ::T
    τ::Vector{T}
    element::Int32 # TODO: deberia llamarlo cell_id o algo asi
end

function Segment(p::Point2D{T}, q::Point2D{T}, element::Int32=Int32(-1)) where {T}
    return Segment(p, q, norm(p - q), Vector{T}(), element)
end

ℓ(segment::Segment) = segment.ℓ

in(x::Point2D, segment::Segment) = point_in_segment(segment.p, segment.q, x)

function point_in_segment(p::Point2D, q::Point2D, x::Point2D)
    lpx = norm(p - x)
    lqx = norm(q - x)
    lpq = norm(p - q)
    return isapprox(lpx + lqx, lpq)
end