
struct Segment{I,O,T}
    xi::I # son Gridap.Point o podrian ser un SVector ya que es semejante y tienen banda de operaciones hechas
    xo::O
    ℓ::T
    τ::T
    # element::E
end

# UID es universal id and AZ es azimuthal_index
# struct Track{UID,AZMI,T,I,O,S}
struct Track{T<:Real,I,O,S}
    ϕ::T
    xi::I
    xo::O
    ℓ::T
    n::S
    segments::Vector{Segment}

    # por ahora comento esto:

    # boundary_i::PhysicalEntity
    # boundary_o::PhysicalEntity

    # next_track_fwd::Track
    # next_track_bwd::Track

    # dir_next_track_fwd::Int
    # dir_next_track_bwd::Int
end
