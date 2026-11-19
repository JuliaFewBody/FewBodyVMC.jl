raw"""
    DiatomicOneElectronWavefunction(ops; α=1.2, κ=6.0, R₀=2.1)

Positive, rotationally invariant trial wavefunction for two nuclei and one
electron,

```math
\psi = e^{-\kappa(R-R_0)^2}\left(e^{-\alpha r_{1e}}+e^{-\alpha r_{2e}}\right).
```

The first three Coulomb weights in `ops` must be ordered as nucleus–nucleus,
first-nucleus–electron, and second-nucleus–electron. This is the order produced
by `Coulomb(charges)` in a shared `Hamiltonian`, or the legacy
`ops += "Coulomb"`, for particle ordering `[nucleus1, nucleus2, electron]`.
"""
struct DiatomicOneElectronWavefunction{T<:Real,V<:AbstractVector{T}}
    w_nn::V
    w_ne1::V
    w_ne2::V
    α::T
    κ::T
    R₀::T
end

function DiatomicOneElectronWavefunction(ops; α::Real = 1.2, κ::Real = 6.0, R₀::Real = 2.1)
    weights = coulomb_weights(ops)
    length(weights) == 3 || throw(
        ArgumentError(
            "a diatomic one-electron system requires exactly three Coulomb pairs",
        ),
    )
    α > 0 || throw(ArgumentError("α must be positive"))
    κ > 0 || throw(ArgumentError("κ must be positive"))
    R₀ > 0 || throw(ArgumentError("R₀ must be positive"))
    return DiatomicOneElectronWavefunction(
        Float64.(weights[1]),
        Float64.(weights[2]),
        Float64.(weights[3]),
        Float64(α),
        Float64(κ),
        Float64(R₀),
    )
end

function _pair_vector(position, w)
    ncoords = length(w)
    length(position) == 3ncoords || throw(
        DimensionMismatch(
            "position must contain three Cartesian components per Jacobi coordinate",
        ),
    )
    X = reshape(position, ncoords, 3)
    return transpose(X) * w
end

_pair_distance(position, w) = norm(_pair_vector(position, w))

function (wavefunction::DiatomicOneElectronWavefunction)(position::AbstractVector)
    R = _pair_distance(position, wavefunction.w_nn)
    r1 = _pair_distance(position, wavefunction.w_ne1)
    r2 = _pair_distance(position, wavefunction.w_ne2)
    nuclear = exp(-wavefunction.κ * (R - wavefunction.R₀)^2)
    return nuclear * (exp(-wavefunction.α * r1) + exp(-wavefunction.α * r2))
end
