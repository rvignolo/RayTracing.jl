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