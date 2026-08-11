```@meta
CurrentModule = FewBodyVMC
```

# FewBodyVMC.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaFewBody.github.io/FewBodyVMC.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaFewBody.github.io/FewBodyVMC.jl/dev/)
[![Build Status](https://github.com/JuliaFewBody/FewBodyVMC.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaFewBody/FewBodyVMC.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaFewBody/FewBodyVMC.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaFewBody/FewBodyVMC.jl)

FewBodyVMC.jl solves quantum few-body bound-state problems using variational Monte Carlo. It combines Hamiltonians from FewBodyHamiltonians.jl with trial wave functions and provides a solve-based interface for energy estimation, sampling diagnostics, and convergence analysis.

## Installation

Run the following command in the Julia REPL or a notebook:

```julia
import Pkg; Pkg.add(url="https://github.com/JuliaFewBody/FewBodyVMC.jl.git")
```

After installation, run the following to load the package and verify it works:

```@repl
import FewBodyVMC; FewBodyVMC.hello()
```

## API Reference

For the generated API index and docstrings, see the [API Reference](api.md).

## Acknowledgments

This package is written in the [Julia programming language](https://julialang.org/), built on an initial project template generated using [PkgFactory.jl](https://github.com/ohno/PkgFactory.jl). This repository is hosted on [GitHub](https://github.com/JuliaFewBody/FewBodyVMC.jl), and continuous integration is run using [GitHub Actions](https://github.com/JuliaFewBody/FewBodyVMC.jl/actions).
