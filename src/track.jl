
struct Point2D{T} <: FieldVector{2,T}
    x::T
    y::T
end

Point2D(x, y) = Point2D(promote(x, y)...)

StaticArrays.similar_type(::Type{<:Point2D}, ::Type{T}, ::Size{(2,)}) where {T} = Point2D{T}

Base.convert(::Type{<:Point2D{T}}, arg::VectorValue{D}) where {D,T} = Point2D{T}(Tuple(arg))

struct Segment{I,O,E,T}
    xi::I
    xo::O
    ℓ::T
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
