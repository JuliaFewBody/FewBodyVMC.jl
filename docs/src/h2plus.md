```@meta
CurrentModule = FewBodyVMC
```

# H₂⁺ without Born–Oppenheimer

This example reproduces the physical setup of the
[FewBodyECG.jl H₂⁺ example](https://juliafewbody.github.io/FewBodyECG.jl/dev/examples/h2plus/)
with variational Monte Carlo. The two protons and the electron are treated as a
direct three-body Coulomb system, including nuclear motion; no
Born–Oppenheimer separation is made.

The nonrelativistic reference energy used by the ECG example is
``E_\mathrm{ref} = -0.597139`` Ha. The dissociation threshold for
``\mathrm{H} + p^+`` is approximately ``-0.5`` Ha, so an energy below that
threshold demonstrates binding.

## Hamiltonian and trial function

The particle order, masses, charges, and operator construction match the ECG
example. FewBodyVMC.jl uses the same mass-scaled Jacobi-coordinate convention,
so the centre-of-mass motion is removed automatically.

The positive, rotationally invariant trial function is

```math
\psi = \exp[-\kappa(R-R_0)^2]
\left(\exp[-\alpha r_{1e}] + \exp[-\alpha r_{2e}]\right),
```

where ``R`` is the proton-proton distance and ``r_{1e}``, ``r_{2e}`` are the
proton-electron distances. It is a compact two-centre ansatz rather than the
systematically enlarged explicitly correlated Gaussian basis used by
FewBodyECG.jl.

```@example h2plus
using FewBodyVMC

mₚ = 1836.15267343
ops = Operators([mₚ, mₚ, 1.0], [+1.0, +1.0, -1.0])
ops += "Kinetic"
ops += "Coulomb"

ψ = DiatomicOneElectronWavefunction(ops; α = 1.2, κ = 6.0, R₀ = 2.1)

# Start with the protons 2 bohr apart. jacobi_coordinates removes the
# centre-of-mass coordinate and returns the package's internal coordinates.
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
    δ = [12.0, 1.6], # proton-proton and electronic Jacobi proposal widths
    r₀ = initial,
    block_size = 50,
)

solution = solve(ops, ψ, method)
solution
```

Compare the estimate with the ECG reference and test the same binding
criterion:

```@example h2plus
h2plus_reference = -0.597139
(
    E₀ = solution.E₀,
    standard_error = solution.standard_error,
    acceptance_rate = solution.acceptance_rate,
    reference = h2plus_reference,
    difference = solution.E₀ - h2plus_reference,
    bound_below_H_plus_p = solution.E₀ < -0.5,
)
```

With the fixed default random-number seed, the estimate should be around
``-0.59`` Ha and the binding check should be `true`. The exact displayed value
is a Monte Carlo estimate: `solution.standard_error` measures sampling
uncertainty using contiguous energy blocks, and `solution.acceptance_rate` is a
diagnostic for the proposal widths.

At this sample size, the difference from ``-0.597139`` Ha contains both
sampling uncertainty and variational bias from the compact trial function.
Increasing `n_steps` reduces the former but does not systematically remove the
latter; doing so requires a more flexible trial function or optimization of
``\alpha``, ``\kappa``, and ``R_0``.
