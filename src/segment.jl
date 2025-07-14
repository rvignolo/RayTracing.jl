"""
    Segment{T<:Real}

Represents a straight-line segment within a single mesh element during ray tracing.

When a [`Track`](@ref) traverses through the computational domain, its intersection with
each mesh element creates a `Segment`. This discretization enables numerical integration
of transport equations over individual mesh cells.

## Fields
- `p::Point2D{T}`: Entry point where the track enters the mesh element
- `q::Point2D{T}`: Exit point where the track exits the mesh element
- `ℓ::T`: Length of the segment within the element
- `τ::Vector{T}`: Storage for transport-related data (e.g., optical thickness)
- `element::Int32`: ID of the mesh element containing this segment

## Usage
Segments are automatically generated during track segmentation and used for
finite element calculations in transport sweeps.

See also: [`Track`](@ref), [`segmentize!`](@ref)
"""
struct Segment{T<:Real}
    p::Point2D{T}
    q::Point2D{T}
    ℓ::T
    τ::Vector{T}
    element::Int32
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