"""
    Segment{I<:Point2D,O<:Point2D,T<:Real}

The intersection of a mesh element or cell and a [`Track`](@ref) yields to a `Segment` with
entry point `xi`, exit point `xo` and length `ℓ`.
"""
struct Segment{I<:Point2D,O<:Point2D,T<:Real}
    xi::I
    xo::O
    ℓ::T
    element::Int32
end

Segment(p::Point2D, q::Point2D, element::Int32=Int32(-1)) = Segment(p, q, norm(p - q), element)

in(x::Point2D, s::Segment) = point_in_segment(s.xi, s.xo, x)

function point_in_segment(p::Point2D, q::Point2D, x::Point2D)
    lpx = norm(p - x)
    lqx = norm(q - x)
    lpq = norm(p - q)
    return isapprox(lpx + lqx, lpq)
end