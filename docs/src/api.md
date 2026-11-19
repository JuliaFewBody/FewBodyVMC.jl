```@meta
CurrentModule = FewBodyVMC
```

# API Reference

The API reference below is generated from docstrings in the [source code](https://github.com/JuliaFewBody/FewBodyVMC.jl/tree/main/src).

## Index

```@index
```

## Docstrings

### Shared Hamiltonian construction

These types are re-exported from FewBodyHamiltonians. VMC accepts shared
Hamiltonians containing atomic-unit `Kinetic` and `Coulomb` terms only; the
remaining construction types are available for use with other backends.

See the [FewBodyHamiltonians API](https://JuliaFewBody.github.io/FewBodyHamiltonians.jl/dev/)
for the shared construction docstrings: `Hamiltonian`, `UnitSystem`, `Kinetic`,
`Coulomb`, `Gaussian`, `ManyBodyGaussian`, `Oscillator`, `Custom`, `Tensor`,
`SpinOrbit`, and the extension supertypes `Operator`, `KineticTerm`, and
`PotentialTerm`.

### VMC

```@autodocs
Modules = [FewBodyVMC]
```
