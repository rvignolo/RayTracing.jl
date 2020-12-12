
# IDEA: use `SArray`s instead? ahi calculo todos los numeros de una en un SVector y ya me
# parece que va a ser mejor, incluso con las quadraturas azimutales (supongo que no me van a
# dar tantos angulos) y las totales
struct PolarQuadrature{N,T<:Real}
    sinθs::Vector{T} # SVector{N,T} como type
    θs::Vector{T}
    ωₚ::Vector{T}
end

npolar(::PolarQuadrature{N}) where {N} = N
npolar2(::PolarQuadrature{N}) where {N} = div(N, 2)

for polar_quadrature in (:TabuchiYamamoto, :GaussLegendre, :EqualWeight, :EqualAngle, :Leonard)
    @eval begin
        $polar_quadrature(N::Int, T::Type{<:Real}=Float64) = $polar_quadrature(Val(N), T)
    end
end

const _TL_NΘ_ = (2, 4, 6)
const _GL_NΘ_ = (2, 3, 4, 6, 8, 10, 12)

# this can be called directly if type stability is needed
function TabuchiYamamoto(n_polar::Val{N}, T::Type{<:Real}=Float64) where {N}
    N in _TL_NΘ_ || throw(DomainError(N, "number of polar angles for TabuchiYamamoto must be in $(string(_TL_NΘ_))."))

    sinθs, θs, ωₚ = ntuple(_ -> Vector{T}(undef, N), 3)

    # Nuclear Engineering Handbook, pg. 1154
    set_tabuchiyamamoto_data!(n_polar, sinθs, ωₚ)

    # compute remaining data from loaded data
    for i in 1:div(N, 2)
        θs[i] = asin(sinθs[i])

        # suplementaries
        j = N - i + 1
        sinθs[j] = sinθs[i]
        θs[j] = π - θs[i]
        ωₚ[j] = ωₚ[i]
    end

    return PolarQuadrature{N,T}(sinθs, θs, ωₚ)
end

function set_tabuchiyamamoto_data!(::Val{2}, sinθs, ωₚ)
    sinθs[1] = 0.798184
    ωₚ[1] = 1 / 2
    return nothing
end

function set_tabuchiyamamoto_data!(::Val{4}, sinθs, ωₚ)
    sinθs[1] = 0.363900
    sinθs[2] = 0.899900
    ωₚ[1] = 0.212854 / 2
    ωₚ[2] = 0.787146 / 2
    return nothing
end

function set_tabuchiyamamoto_data!(::Val{6}, sinθs, ωₚ)
    sinθs[1] = 0.166648
    sinθs[2] = 0.537707
    sinθs[3] = 0.932954
    ωₚ[1] = 0.046233 / 2
    ωₚ[2] = 0.283619 / 2
    ωₚ[3] = 0.670148 / 2
    return nothing
end

function GaussLegendre(n_polar::Val{N}, T::Type{<:Real}=Float64) where {N}
    N in _GL_NΘ_ || throw(DomainError(N, "number of polar angles for GaussLegendre must be in $(string(_GL_NΘ_))."))

    sinθs, θs, ωₚ = ntuple(_ -> Vector{T}(undef, N), 3)

    # Nuclear Engineering Handbook, pg. 1153
    set_gausslegendre_data!(n_polar, θs, ωₚ)

    # compute remaining data from loaded data
    for i in 1:div(N, 2)
        sinθs[i] = sin(θs[i])

        # suplementaries
        j = N - i + 1
        sinθs[j] = sinθs[i]
        θs[j] = π - θs[i]
        ωₚ[j] = ωₚ[i]
    end

    return PolarQuadrature{N,T}(sinθs, θs, ωₚ)
end

function set_gausslegendre_data!(::Val{2}, θs, ωₚ)
    θs[1] = acos(0.5773502691)
    ωₚ[1] = 1 / 2
    return nothing
end

function set_gausslegendre_data!(::Val{4}, θs, ωₚ)
    θs[1] = acos(0.3399810435)
    θs[2] = acos(0.8611363115)
    ωₚ[1] = 0.6521451549 / 2
    ωₚ[2] = 0.3478548451 / 2
    return nothing
end

function set_gausslegendre_data!(::Val{6}, θs, ωₚ)
    θs[1] = acos(0.2386191860)
    θs[2] = acos(0.6612093864)
    θs[3] = acos(0.9324695142)
    ωₚ[1] = 0.4679139346 / 2
    ωₚ[2] = 0.3607615730 / 2
    ωₚ[3] = 0.1713244924 / 2
    return nothing
end

function set_gausslegendre_data!(::Val{8}, θs, ωₚ)
    θs[1] = acos(0.1834346424)
    θs[2] = acos(0.5255324099)
    θs[3] = acos(0.7966664774)
    θs[4] = acos(0.9602898564)
    ωₚ[1] = 0.3626837834 / 2
    ωₚ[2] = 0.3137066459 / 2
    ωₚ[3] = 0.2223810344 / 2
    ωₚ[4] = 0.1012285363 / 2
    return nothing
end

function set_gausslegendre_data!(::Val{10}, θs, ωₚ)
    θs[1] = acos(0.1488743387)
    θs[2] = acos(0.4333953941)
    θs[3] = acos(0.6794095682)
    θs[4] = acos(0.8650633666)
    θs[5] = acos(0.9739065285)
    ωₚ[1] = 0.2955242247 / 2
    ωₚ[2] = 0.2692667193 / 2
    ωₚ[3] = 0.2190863625 / 2
    ωₚ[4] = 0.1494513492 / 2
    ωₚ[5] = 0.0666713443 / 2
    return nothing
end

function set_gausslegendre_data!(::Val{12}, θs, ωₚ)
    θs[1] = acos(0.1252334085)
    θs[2] = acos(0.3678314989)
    θs[3] = acos(0.5873179542)
    θs[4] = acos(0.7699026741)
    θs[5] = acos(0.9041172563)
    θs[6] = acos(0.9815606342)
    ωₚ[1] = 0.2491470458 / 2
    ωₚ[2] = 0.2334925365 / 2
    ωₚ[3] = 0.2031674267 / 2
    ωₚ[4] = 0.1600783286 / 2
    ωₚ[5] = 0.1069393260 / 2
    ωₚ[6] = 0.0471753364 / 2
    return nothing
end

function Leonard(n_polar::Val{N}, T::Type{<:Real}=Float64) where {N}
    N in _TL_NΘ_ || throw(DomainError(N, "number of polar angles for Leonard must be in $(string(_TL_NΘ_))."))

    sinθs, θs, ωₚ = ntuple(_ -> Vector{T}(undef, N), 3)

    # TODO: decide which source we should use since Tabuchi paper and OpenMoC have differ.
    # On the other hand, the Nuclear Engineering Handbook has Herbert optimized values.
    # The following values are from Tabuchi paper
    set_leonard_data!(n_polar, θs, ωₚ)

    # compute remaining data from loaded data
    for i in 1:div(N, 2)
        θs[i] = asin(sinθs[i])

        # suplementaries
        j = N - i + 1
        sinθs[j] = sinθs[i]
        θs[j] = π - θs[i]
        ωₚ[j] = ωₚ[i]
    end

    return PolarQuadrature{N,T}(sinθs, θs, ωₚ)
end

function set_leonard_data!(::Val{2}, sinθs, ωₚ)
    sinθs[1] = 0.752244
    ωₚ[1] = 1 / 2
    return nothing
end

function set_leonard_data!(::Val{4}, sinθs, ωₚ)
    sinθs[1] = 0.273658
    sinθs[2] = 0.865714
    ωₚ[1] = 0.139473 / 2
    ωₚ[2] = 0.860527 / 2
    return nothing
end

function set_leonard_data!(::Val{6}, sinθs, ωₚ)
    sinθs[1] = 0.103840
    sinθs[2] = 0.430723
    sinθs[3] = 0.905435
    ωₚ[1] = 0.020530 / 2
    ωₚ[2] = 0.219161 / 2
    ωₚ[3] = 0.760309 / 2
    return nothing
end

struct AzimuthalQuadrature{T}
    n_azim::Int    # the number of azimuthal angles in (0, 2π)
    n_azim_2::Int  # the number of azimuthal angles in (0, π)
    n_azim_4::Int  # the number of azimuthal angles in (0, π/2)

    δ::T # azimuthal spacing (provided as input)
    δs::Vector{T} # effective azimuthal spacings

    ϕs::Vector{T}
    ωₐ::Vector{T}
end

function AzimuthalQuadrature(n_azim, δ::T) where {T<:Real}
    n_azim_2 = div(n_azim, 2)
    n_azim_4 = div(n_azim, 4)

    δs, ϕs, ωₐ = ntuple(_ -> Vector{T}(undef, n_azim_2), 3)

    return AzimuthalQuadrature(n_azim, n_azim_2, n_azim_4, δ, δs, ϕs, ωₐ)
end

function init_weights!(aq::AzimuthalQuadrature)
    @unpack ϕs, ωₐ, n_azim_4 = aq
    for i in right_dir(aq)
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

right_dir(aq::AzimuthalQuadrature) = 1:aq.n_azim_4
left_dir(aq::AzimuthalQuadrature) = (aq.n_azim_4+1):aq.n_azim_2
both_dir(aq::AzimuthalQuadrature) = 1:aq.n_azim_2

points_right(aq::AzimuthalQuadrature, i) = i <= aq.n_azim_4
points_left(aq::AzimuthalQuadrature, i) = !points_right(aq, i)

suplementary_idx(aq::AzimuthalQuadrature, i) = aq.n_azim_2 - i + 1

struct Quadrature{A<:AzimuthalQuadrature,P<:PolarQuadrature,T<:Real}
    azimuthal::A
    polar::P
    ω::Matrix{T}
end

function Quadrature(azimuthal::AzimuthalQuadrature{T}, polar::PolarQuadrature{N,T}, i) where {N,T}
    @unpack n_azim_2, δs, ωₐ = azimuthal
    @unpack sinθs, ωₚ = polar
    n_polar_2 = npolar2(polar)

    # IDEA: I think we can use n_azim_4 because the matrix shows repeated values
    ω = Matrix{T}(undef, n_azim_2, n_polar_2) # we have polar symmetry

    for i in both_dir(azimuthal), j in 1:n_polar_2
        ω[i, j] = 4 * π * ωₐ[i] * ωₚ[j] * δs[i] * sinθs[j]
    end
    # ω .*= 4π

    return Quadrature(azimuthal, polar, ω)
end