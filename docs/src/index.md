```@meta
CurrentModule = FewBodyVMC
```

# FewBodyVMC.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaFewBody.github.io/FewBodyVMC.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaFewBody.github.io/FewBodyVMC.jl/dev/)
[![Build Status](https://github.com/JuliaFewBody/FewBodyVMC.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaFewBody/FewBodyVMC.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaFewBody/FewBodyVMC.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaFewBody/FewBodyVMC.jl)

FewBodyVMC estimates nonrelativistic few-body ground-state energies with
variational Monte Carlo. It uses the same mass-scaled Jacobi coordinates and
particle-indexed Hamiltonians from
[FewBodyHamiltonians.jl](https://github.com/JuliaFewBody/FewBodyHamiltonians.jl)
as [FewBodyECG.jl](https://github.com/JuliaFewBody/FewBodyECG.jl).

Start with the [HD⁺ example](example.md), which treats the proton, deuteron,
and electron as a full three-body Coulomb system without a Born–Oppenheimer
approximation. The [H₂⁺ example](h2plus.md) reproduces the physical system and
binding check from the corresponding FewBodyECG.jl example with VMC.

## Installation

Run the following command in the Julia REPL or a notebook:

```julia
import Pkg; Pkg.add(url="https://github.com/JuliaFewBody/FewBodyVMC.jl.git")
```

After installation, load the package and construct the HD⁺ Hamiltonian:

```julia
using FewBodyVMC

H = Kinetic([1836.15267343, 3670.48296788, 1.0]) + Coulomb([1.0, 1.0, -1.0])
```

`Kinetic` defines particle masses and removes center-of-mass motion;
`Coulomb(charges)` adds all pairs in particle order. These constructors and
`Hamiltonian` are re-exported from FewBodyHamiltonians without wrapping their
types. Pass the same `H` to trial-wavefunction construction, coordinate helpers,
`local_energy`, and `solve`.

VMC currently supports atomic-unit `Kinetic` and `Coulomb` terms only. Other
shared potentials (including `Custom`) and non-atomic units fail before
wavefunction evaluation or sampling. The legacy `Operators(masses, charges)`
builder and string additions remain available.

## API Reference

For the generated API index and docstrings, see the [API Reference](api.md).

## Acknowledgments

This package is written in the [Julia programming language](https://julialang.org/), built on an initial project template generated using [PkgFactory.jl](https://github.com/ohno/PkgFactory.jl). This repository is hosted on [GitHub](https://github.com/JuliaFewBody/FewBodyVMC.jl), and continuous integration is run using [GitHub Actions](https://github.com/JuliaFewBody/FewBodyVMC.jl/actions).
