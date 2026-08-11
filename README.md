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

After installation, load the package and verify it works:

```julia
julia> import FewBodyVMC; FewBodyVMC.hello()
"Hello, World!"
```

## Documentation

- Home: https://JuliaFewBody.github.io/FewBodyVMC.jl
- API Reference: https://JuliaFewBody.github.io/FewBodyVMC.jl/dev/api
