
mutable struct PolarQuadrature{T<:Real}
    n_polar::Int
    n_polar_2::Int

    sinθs::Vector{T}
    θs::Vector{T}
    ωₚ::Vector{T}
end

struct Quadrature{T<:Real,Q<:PolarQuadrature}
    n_azim::Int    # el numero de angulos azimutales en (0, 2π)
    n_azim_2::Int  # el numero de angulos azimutales en (0, π)
    n_azim_4::Int  # el numero de angulos azimutales en (0, π/2)

    δ::T # azimuthal spacing
    δs::Vector{T} # effective azimuthal spacings

    ϕs::Vector{T}
    ωₐ::Vector{T}

    polarquad::Q

    ω::Matrix{T}   # total weight, creo que basta con que sea matrix porque es rectangular
end