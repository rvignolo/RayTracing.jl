
struct Point2D{T} <: FieldVector{2,T}
    x::T
    y::T
end

StaticArrays.similar_type(::Type{<:Point2D}, ::Type{T}, s::Size{(2,)}) where {T} = Point2D{T}

struct Segment{I,O,E,T}
    xi::I
    xo::O
    ℓ::T
    # τ::T
    element::E
end

struct Track{T<:Real,I,O,S}
    ϕ::T
    xi::I
    xo::O
    ℓ::T
    n::S
    segments::Vector{Segment}

    # we will eventually need the following fields fields:

    # boundary_i::PhysicalEntity
    # boundary_o::PhysicalEntity

    # next_track_fwd::Track
    # next_track_bwd::Track

    # dir_next_track_fwd::Int
    # dir_next_track_bwd::Int
end
