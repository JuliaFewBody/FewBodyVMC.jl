raw"""
    DiatomicOneElectronWavefunction(ops; α=1.2, κ=6.0, R₀=2.1)

Positive, rotationally invariant trial wavefunction for two nuclei and one
electron,

```math
\psi = e^{-\kappa(R-R_0)^2}\left(e^{-\alpha r_{1e}}+e^{-\alpha r_{2e}}\right).
```

For a shared `Hamiltonian`, particles 1 and 2 are the nuclei and particle 3 is
the electron. Include all three Coulomb pairs in any order. For Jacobi-frame
operators, supply the three weights in the order (1,2), (1,3), (2,3), as
produced by `ops += "Coulomb"`.
"""
struct DiatomicOneElectronWavefunction{T<:Real,V<:AbstractVector{T}}
    w_nn::V
    w_ne1::V
    w_ne2::V
    α::T
    κ::T
    R₀::T
end

function _diatomic_weights(ops)
    weights = coulomb_weights(ops)
    length(weights) == 3 || throw(
        ArgumentError(
            "a diatomic one-electron system requires exactly three Coulomb pairs",
        ),
    )
    return weights
end

function _diatomic_weights(H::Hamiltonian)
    masses = _vmc_kinetic(H).masses
    length(masses) == 3 ||
        throw(ArgumentError("the diatomic trial requires three particles"))
    pairs = Set(minmax(term.i, term.j) for term in H if term isa Coulomb)
    length(pairs) == 3 ||
        throw(ArgumentError("the diatomic trial requires all three Coulomb pairs"))
    _, U = jacobi_transform(masses)
    return (_pair_weight(U, 1, 2), _pair_weight(U, 1, 3), _pair_weight(U, 2, 3))
end

function DiatomicOneElectronWavefunction(ops; α::Real = 1.2, κ::Real = 6.0, R₀::Real = 2.1)
    weights = _diatomic_weights(ops)
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
