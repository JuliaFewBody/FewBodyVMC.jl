module FewBodyVMC

import FewBodyHamiltonians
using FewBodyHamiltonians: Operator, KineticTerm, PotentialTerm, Hamiltonian,
    UnitSystem, Kinetic, Coulomb, Gaussian, ManyBodyGaussian, Oscillator, Custom,
    Tensor, SpinOrbit
using ForwardDiff
using LinearAlgebra
using Random

export Operator,
    KineticTerm,
    PotentialTerm,
    Hamiltonian,
    UnitSystem,
    Kinetic,
    Coulomb,
    Gaussian,
    ManyBodyGaussian,
    Oscillator,
    Custom,
    Tensor,
    SpinOrbit,
    Operators,
    KineticOperator,
    CoulombOperator,
    VariationalMonteCarlo,
    VMCSolution,
    DiatomicOneElectronWavefunction,
    solve,
    local_energy,
    energies,
    jacobi_transform,
    jacobi_coordinates,
    particle_coordinates,
    coulomb_weights,
    Λ

include("coordinates.jl")
include("operators.jl")
include("wavefunctions.jl")
include("solution.jl")
include("vmc.jl")

end
