```@meta
CurrentModule = FewBodyVMC
```

# Example

## HD⁺ ground state

This example treats HD⁺ as the nonrelativistic three-body Coulomb system
``(p^+, d^+, e^-)`` in atomic units. The centre-of-mass motion is removed with
the same mass-scaled Jacobi transformation used by FewBodyECG.jl.

The positive, rotationally invariant trial function is

```math
\psi = \exp[-\kappa(R-R_0)^2]
\left(\exp[-\alpha r_{pe}] + \exp[-\alpha r_{de}]\right).
```

It represents a nuclear vibrational envelope multiplying the simplest
two-centre molecular orbital. The parameters below are useful starting values,
not a high-precision optimization.

```@example hdplus
using FewBodyVMC

# Particle order is proton, deuteron, electron. Masses and energies are in
# electron masses and hartree, respectively.
mₚ = 1836.15267343
m_d = 3670.48296788
ops = Operators([mₚ, m_d, 1.0], [+1.0, +1.0, -1.0])
ops += "Kinetic"
ops += "Coulomb"

ψ = DiatomicOneElectronWavefunction(ops; α = 1.2, κ = 6.0, R₀ = 2.1)

# Start near a 2-bohr internuclear separation. Conversion to the internal
# coordinates is explicit, so particle and Jacobi conventions cannot be mixed.
particle_positions = [
    -1.0  0.0  0.0
     1.0  0.0  0.0
     0.0  0.2  0.0
]
initial = vec(jacobi_coordinates(ops, particle_positions))

method = VariationalMonteCarlo(
    n_steps = 2_500,
    burn_in = 750,
    thinning = 2,
    n_walkers = 4,
    δ = [12.0, 1.6], # nuclear and electronic Jacobi proposal widths
    r₀ = initial,
    block_size = 50,
)

solution = solve(ops, ψ, method)
solution
```

The result should fluctuate near ``-0.58`` Ha and below the ``\mathrm{H}+d^+``
dissociation threshold (approximately ``-0.5`` Ha), demonstrating a bound
ground-state variational estimate. The accurate nonrelativistic value is about
``-0.59789797`` Ha [1]; the difference is primarily variational bias from the
compact trial function. Increasing only `n_steps` reduces statistical error,
not that bias. A richer trial function or parameter optimization is the next
accuracy step.

The two entries in `δ` matter: mass scaling makes the nuclear Jacobi coordinate
much broader than the electronic coordinate. `standard_error` is computed from
contiguous energy blocks, while `acceptance_rate` is useful for retuning the
proposal widths.

1. L. Hilico, N. Billy, B. Grémaud, and D. Delande,
   [*Ab initio calculation of the J = 0 and J = 1 states of the H₂⁺, D₂⁺ and
   HD⁺ molecular ions*](https://citeseerx.ist.psu.edu/document?doi=ec487a35ae167f493613639759170ac1c95d1271),
   *European Physical Journal D* **12**, 449–466 (2000).
