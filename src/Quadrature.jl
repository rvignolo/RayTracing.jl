
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

# this can be called if type stability is needed
function TabuchiYamamoto(n_polar::Val{N}, T::Type{<:Real}=Float64) where {N}
    @assert N in (2, 4, 6) "the number of polar angles for TabuchiYamamoto must be in (2, 4, 6)"

    sinθs, θs, ωₚ = ntuple(_ -> Vector{T}(undef, N), 3)

    # Nuclear Engineering Handbook, pg. 1154
    set_tabuchiyamamoto_data!(n_polar, sinθs, ωₚ)

    # a partir de la data, calculamos los datos que faltan
    for i in 1:div(N, 2)
        j = N - i + 1
        sinθs[j] = sinθs[i]
        ωₚ[j] = ωₚ[i]
        θs[i] = asin.(sinθs[i])
        θs[j] = π - θs[i]
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

function GaussLegendre(n_polar::Val{N}, T::Type{<:Real}=Float64) where {N<:Int}
    # filter(iseven, 1:12)
    @assert N in (2, 4, 6, 8, 10, 12) "the number of polar angles for TabuchiYamamoto must be in (2, 4, 6, 8, 10, 12)"

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