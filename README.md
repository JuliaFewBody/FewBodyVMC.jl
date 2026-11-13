# FewBodyVMC.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaFewBody.github.io/FewBodyVMC.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaFewBody.github.io/FewBodyVMC.jl/dev/)
[![Build Status](https://github.com/JuliaFewBody/FewBodyVMC.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaFewBody/FewBodyVMC.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaFewBody/FewBodyVMC.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaFewBody/FewBodyVMC.jl)

FewBodyVMC.jl solves quantum few-body bound-state problems using variational Monte Carlo. It combines Hamiltonians from FewBodyHamiltonians.jl with trial wave functions and provides a solve-based interface for energy estimation, sampling diagnostics, and convergence analysis.

## Quick Start

Run the following command in the Julia REPL or a notebook:

```julia
import Pkg; Pkg.add(url="https://github.com/JuliaFewBody/FewBodyVMC.jl.git")
```

After installation, build a three-body Hamiltonian and trial wavefunction:

```julia
using FewBodyVMC

H = Kinetic([1836.15267343, 3670.48296788, 1.0]) + Coulomb([1.0, 1.0, -1.0])
ψ = DiatomicOneElectronWavefunction(H)
```

`Hamiltonian`, `Kinetic`, and `Coulomb` are the exact constructors from
[FewBodyHamiltonians.jl](https://github.com/JuliaFewBody/FewBodyHamiltonians.jl),
also used by FewBodyECG. VMC currently accepts atomic-unit kinetic and Coulomb
terms; unsupported potentials and unit systems are rejected before sampling.
The existing `Operators(masses, charges)` builder with `ops += "Kinetic"` and
`ops += "Coulomb"` remains supported.

## Documentation

- Home: https://JuliaFewBody.github.io/FewBodyVMC.jl
- HD⁺ Example: https://JuliaFewBody.github.io/FewBodyVMC.jl/dev/example
- H₂⁺ Example: https://JuliaFewBody.github.io/FewBodyVMC.jl/dev/h2plus
- API Reference: https://JuliaFewBody.github.io/FewBodyVMC.jl/dev/api

## Citation

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
