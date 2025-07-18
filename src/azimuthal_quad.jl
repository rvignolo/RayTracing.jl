"""
    AzimuthalQuadrature{N,N2,N4,T<:Real}

Represents the angular discretization for ray tracing in neutron transport calculations.

This structure holds all information related to the azimuthal quadrature, including the discretized
azimuthal angles, their associated weights for numerical integration, and the spacing parameters
used in the ray tracing algorithm.

## Type Parameters
- `N`: Total number of azimuthal angles (must be a multiple of 4)
- `N2`: Number of angles in the half-plane (0, π), equal to `N÷2`
- `N4`: Number of angles in the first quadrant (0, π/2), equal to `N÷4`
- `T<:Real`: Numeric type for all calculations (Float64, Float32, etc.)

## Fields
- `δ::T`: Nominal azimuthal spacing used to generate the quadrature
- `δs::Vector{T}`: Effective azimuthal spacings for each angle (length N2)
- `ϕs::Vector{T}`: Azimuthal angles in radians (length N2, spanning 0 to π)
- `ωₐ::Vector{T}`: Angular weights for numerical integration (length N2)

## Mathematical Structure

The azimuthal quadrature organizes angles into a symmetric structure:

- **Full domain**: [0, 2π) with N angles total
- **Half-plane**: [0, π) with N2 angles (unique directions)
- **First quadrant**: [0, π/2) with N4 angles
- **Second quadrant**: [π/2, π) with N4 angles

Angles in (π, 2π) are avoided because of the symmetry of the problem. This helps reducing the
computational cost.

## Usage

```julia
# Create quadrature with 8 angles and spacing 0.02
aq = AzimuthalQuadrature(8, 0.02)

# Access quadrature properties
n_total = n_azim_total(aq)      # 8
n_half = n_azim_half(aq)        # 4
n_quad = n_azim_quad(aq)        # 2

# Get angle ranges
quad1 = azimuthal_quadrant_1(aq)  # 1:2
quad2 = azimuthal_quadrant_2(aq)  # 3:4

# Check direction
is_rightward_direction(aq, 1)     # true
is_leftward_direction(aq, 3)      # true

# Find supplementary angle
supp_idx = supplementary_azimuthal_idx(aq, 1)  # 4
```

## Integration Weights

The weights `ωₐ` are computed to ensure proper numerical integration:

```julia
# Sum of weights should integrate to unity
sum(aq.ωₐ) ≈ 0.5  # Over half-plane (0, π)
```

## Notes

- The structure is designed for 2D ray tracing in neutron transport
- Angles are automatically organized to exploit symmetry
- Weights are normalized for proper integration over the angular domain
- The quadrature supports various numeric types for precision control

See also: [`init_weights!`](@ref), [`n_azim_total`](@ref), [`azimuthal_quadrant_1`](@ref)
"""
struct AzimuthalQuadrature{N,N2,N4,T<:Real}
    δ::T
    δs::Vector{T}
    ϕs::Vector{T}
    ωₐ::Vector{T}
end

"""
    n_azim_total(aq::AzimuthalQuadrature)

Returns the total number of azimuthal angles in the quadrature.

This is the full number of azimuthal directions used in the discretization, spanning the complete
angular domain [0, 2π).

## Arguments
- `aq`: Azimuthal quadrature structure

## Returns
- Total number of azimuthal angles (N)

## Example
```julia
aq = AzimuthalQuadrature(8, 0.1)
n_azim_total(aq)  # Returns 8
```
"""
n_azim_total(::AzimuthalQuadrature{N}) where {N} = N

"""
    n_azim_half(aq::AzimuthalQuadrature)

Returns the number of azimuthal angles in the half-plane (0, π).

This represents the number of unique azimuthal directions needed for 2D ray tracing, since angles in
(π, 2π) are supplementary to those in (0, π).

## Arguments
- `aq`: Azimuthal quadrature structure

## Returns
- Number of azimuthal angles in half-plane (N2 = N/2)

## Example
```julia
aq = AzimuthalQuadrature(8, 0.1)
n_azim_half(aq)  # Returns 4
```
"""
n_azim_half(::AzimuthalQuadrature{N,N2}) where {N,N2} = N2

"""
    n_azim_quad(aq::AzimuthalQuadrature)

Returns the number of azimuthal angles in the first quadrant (0, π/2).

This represents the number of unique directions in the first quadrant, which is used for computing
track origins and directions in the ray tracing algorithm.

## Arguments
- `aq`: Azimuthal quadrature structure

## Returns
- Number of azimuthal angles in first quadrant (N4 = N/4)

## Example
```julia
aq = AzimuthalQuadrature(8, 0.1)
n_azim_quad(aq)  # Returns 2
```
"""
n_azim_quad(::AzimuthalQuadrature{N,N2,N4}) where {N,N2,N4} = N4

"""
    azimuthal_quadrant_1(aq::AzimuthalQuadrature)

Returns the indices for the first quadrant of azimuthal angles (0, π/2).

In azimuthal quadrature, angles are divided into quadrants. The first quadrant contains angles from
0 to π/2 radians, corresponding to indices 1 to N4.
"""
azimuthal_quadrant_1(::AzimuthalQuadrature{N,N2,N4}) where {N,N2,N4} = 1:N4

"""
    azimuthal_quadrant_2(aq::AzimuthalQuadrature)

Returns the indices for the second quadrant of azimuthal angles (π/2, π).

The second quadrant contains angles from π/2 to π radians, corresponding to indices N4+1 to N2.
"""
azimuthal_quadrant_2(::AzimuthalQuadrature{N,N2,N4}) where {N,N2,N4} = (N4+1):N2

"""
    azimuthal_half_plane(aq::AzimuthalQuadrature)

Returns the indices for all azimuthal angles in the half-plane (0, π).

This includes both quadrants and corresponds to indices 1 to N2.
"""
azimuthal_half_plane(::AzimuthalQuadrature{N,N2}) where {N,N2} = 1:N2

"""
    is_rightward_direction(aq::AzimuthalQuadrature, i)

Checks if the azimuthal angle at index `i` points in the rightward direction.

An angle points rightward if it has a positive x-component, which corresponds to angles in the first
quadrant (indices 1 to N4).
"""
is_rightward_direction(::AzimuthalQuadrature{N,N2,N4}, i) where {N,N2,N4} = i ≤ N4

"""
    is_leftward_direction(aq::AzimuthalQuadrature, i)

Checks if the azimuthal angle at index `i` points in the leftward direction.

An angle points leftward if it has a negative x-component, which corresponds to angles in the second
quadrant (indices N4+1 to N2).
"""
is_leftward_direction(aq::AzimuthalQuadrature, i) = !is_rightward_direction(aq, i)

"""
    supplementary_azimuthal_idx(aq::AzimuthalQuadrature, i)

Returns the index of the supplementary azimuthal angle for index `i`.

In azimuthal quadrature, angles are organized in pairs where each angle has a supplementary angle
that differs by π radians (180°). This function computes the index of the supplementary angle using
the relationship: `supplementary_idx = N2 - i + 1`, where `N2` is the number of azimuthal angles in
the half-plane (0, π).

## Arguments
- `aq`: Azimuthal quadrature structure
- `i`: Index of the azimuthal angle

## Returns
- Index of the supplementary azimuthal angle

## Example
```julia
# For N=8 azimuthal angles, N2=4
# supplementary_azimuthal_idx(aq, 1) == 4  # angles 1 and 4 are supplementary
# supplementary_azimuthal_idx(aq, 2) == 3  # angles 2 and 3 are supplementary
```
"""
supplementary_azimuthal_idx(::AzimuthalQuadrature{N,N2}, i) where {N,N2} = N2 - i + 1

"""
    AzimuthalQuadrature(n_azim::Int, δ::Real)

Construct an azimuthal quadrature structure for ray tracing with `N` azimuthal angles and spacing
`δ`.

## Arguments
- `N`: Number of azimuthal angles (must be positive and a multiple of 4)
- `δ`: Azimuthal spacing (must be positive)

## Returns
- `AzimuthalQuadrature{N,T,N2,N4}` where `N2 = N÷2` (half-plane angles) and `N4 = N÷4` (quadrant
  angles)

## Validation
- Throws `DomainError` if `N ≤ 0`
- Throws `DomainError` if `N` is not a multiple of 4
- Throws `DomainError` if `δ ≤ 0`

## Notes
The constructor initializes uninitialized vectors for effective spacings (`δs`), azimuthal angles
(`ϕs`), and weights (`ωₐ`), each of length `N2`. These arrays are populated later during the tracing
process.

## Example
```julia
aq = AzimuthalQuadrature(8, 0.02)  # 8 angles, spacing 0.02
```
"""
AzimuthalQuadrature(n_azim::Int, δ::Real) = AzimuthalQuadrature(Val(n_azim), δ)

function AzimuthalQuadrature(::Val{N}, δ::T) where {N,T<:Real}
    N > 0 || throw(DomainError(N, "number of azimuthal angles must be positive."))
    iszero(rem(N, 4)) || throw(DomainError(N, "number of azimuthal angles must be a " *
                                              "multiple of 4."))
    δ > 0 || throw(DomainError(δ, "azimuthal spacing must be positive."))

    N2 = div(N, 2)
    N4 = div(N, 4)

    δs, ϕs, ωₐ = ntuple(_ -> Vector{T}(undef, N2), 3)

    return AzimuthalQuadrature{N,N2,N4,T}(δ, δs, ϕs, ωₐ)
end

"""
    init_weights!(aq::AzimuthalQuadrature)

Initialize the azimuthal quadrature weights `ωₐ` based on the computed azimuthal angles `ϕs`.

This function computes the angular weights for numerical integration over the azimuthal domain. The
weights are calculated using a finite difference approximation of the angular spacing between
adjacent azimuthal angles.

## Weight Calculation

For each azimuthal angle index `i` in the first quadrant:

- **First angle (i=1)**: `ωₐ[i] = (ϕs[i+1] - ϕs[i]) / (4π)`
- **Last angle (i=N4)**: `ωₐ[i] = (π - ϕs[i] - ϕs[i-1]) / (4π)`
- **Middle angles**: `ωₐ[i] = (ϕs[i+1] - ϕs[i-1]) / (4π)`

The weights are then copied to the supplementary angles using the relationship `ωₐ[j] = ωₐ[i]` where
`j = supplementary_azimuthal_idx(aq, i)`.

## Arguments
- `aq`: Azimuthal quadrature structure to initialize

## Returns
- `nothing` (modifies `aq.ωₐ` in-place)

## Mathematical Background

The weights represent the angular measure associated with each azimuthal direction, normalized by
`4π` to ensure proper integration over the full solid angle. This normalization accounts for the
fact that we're working in 2D (azimuthal angles only) but the weights should integrate to unity over
the full angular domain.

## Example
```julia
aq = AzimuthalQuadrature(8, 0.1)
# ... compute ϕs ...
init_weights!(aq)
# Now aq.ωₐ contains the normalized angular weights
```
"""
function init_weights!(aq::AzimuthalQuadrature)
    @unpack ϕs, ωₐ = aq
    n_azim_4 = n_azim_quad(aq)
    inv_4π = inv(4π)

    # Compute weights for the first quadrant
    for i in 1:n_azim_4
        # Calculate weight based on position in the quadrant
        if isone(i)
            # First angle: weight from current to next angle
            ωₐ[i] = ϕs[i+1] - ϕs[i]
        elseif isequal(i, n_azim_4)
            # Last angle: weight from previous angle to π
            ωₐ[i] = π - ϕs[i] - ϕs[i-1]
        else
            # Middle angles: weight spans from previous to next angle
            ωₐ[i] = ϕs[i+1] - ϕs[i-1]
        end

        # Normalize by 4π for proper integration
        ωₐ[i] *= inv_4π

        # Copy weight to supplementary angle (symmetry)
        j = supplementary_azimuthal_idx(aq, i)
        ωₐ[j] = ωₐ[i]
    end

    return nothing
end
