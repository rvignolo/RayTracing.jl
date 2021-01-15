"""
    AzimuthalQuadrature{N,T<:Real}

Holds information of the azimuthal quadrature, such as the azimuthal angles `ϕs`, the
azimuthal spacings `δ` as well as effective azimuthal spacings `δs` and the azimuthal
weights `ωₐ`.
"""
struct AzimuthalQuadrature{N,T<:Real,N2,N4}
    δ::T
    δs::Vector{T}
    ϕs::Vector{T}
    ωₐ::Vector{T}
end

nazim(::AzimuthalQuadrature{N}) where {N} = N
nazim2(::AzimuthalQuadrature{N,T,N2}) where {N,T,N2} = N2
nazim4(::AzimuthalQuadrature{N,T,N2,N4}) where {N,T,N2,N4} = N4

AzimuthalQuadrature(n_azim::Int, δ::Real) = AzimuthalQuadrature(Val(n_azim), δ)

function AzimuthalQuadrature(::Val{N}, δ::T) where {N,T<:Real}
    N > 0 || throw(DomainError(N, "number of azimuthal angles must be positive."))
    iszero(rem(N, 4)) || throw(DomainError(N, "number of azimuthal angles must be a " *
                                              "multiple of 4."))
    δ > 0 || throw(DomainError(δ, "azimuthal spacing must be positive."))

    N2 = div(N, 2)
    N4 = div(N, 4)

    δs, ϕs, ωₐ = ntuple(_ -> Vector{T}(undef, N2), 3)

    return AzimuthalQuadrature{N,T,N2,N4}(δ, δs, ϕs, ωₐ)
end

function init_weights!(aq::AzimuthalQuadrature)
    @unpack ϕs, ωₐ = aq
    n_azim_4 = nazim4(aq)

    for i in 1:n_azim_4
        if isone(i)
            ωₐ[i] = ϕs[i+1] - ϕs[i]
        elseif isequal(i, n_azim_4)
            ωₐ[i] = π - ϕs[i] - ϕs[i-1]
        else
            ωₐ[i] = ϕs[i+1] - ϕs[i-1]
        end
        ωₐ[i] /= 4π

        j = suplementary_idx(aq, i)
        ωₐ[j] = ωₐ[i]
    end
    return nothing
end

# TODO: use north-west y cosas asi?
right_dir(::AzimuthalQuadrature{N,T,N2,N4}) where {N,T,N2,N4} = 1:N4
left_dir(::AzimuthalQuadrature{N,T,N2,N4}) where {N,T,N2,N4} = (N4+1):N2
both_dir(::AzimuthalQuadrature{N,T,N2}) where {N,T,N2} = 1:N2

points_right(::AzimuthalQuadrature{N,T,N2,N4}, i) where {N,T,N2,N4} = i ≤ N4
points_left(aq::AzimuthalQuadrature, i) = !points_right(aq, i)

suplementary_idx(::AzimuthalQuadrature{N,T,N2}, i) where {N,T,N2} = N2 - i + 1