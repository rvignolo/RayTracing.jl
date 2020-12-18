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

point_in_segment(s::Segment, x::Point2D) = point_in_segment(s.xi, s.xo, x)

function point_in_segment(xa::Point2D, xb::Point2D, x::Point2D)
    la = norm(xa - x)
    lb = norm(xb - x)
    ab = norm(xa - xb)
    return isapprox(la + lb, ab)
end