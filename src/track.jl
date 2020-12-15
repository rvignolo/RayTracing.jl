
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

# maybe we need this for the MoC solver in the future
# struct Track{ID,AZ,T,I,O,S} ... end # ID is universal id and AZ es azimuthal_index
# id(::Track{ID}) where {ID} = ID
# azimuthal_idx(::Track{ID,AZ}) where {ID,AZ} = AZ
# Track{ID,AZ}(ϕ, xi, xo, ℓ, n, segments) = ...
