import FewBodyHamiltonians

@testset "shared Hamiltonian construction" begin
    @test isempty(Test.detect_ambiguities(FewBodyVMC, FewBodyHamiltonians; recursive = true))
    for name in (
            :Operator, :KineticTerm, :PotentialTerm, :Hamiltonian, :UnitSystem,
            :Kinetic, :Coulomb, :Gaussian, :ManyBodyGaussian, :Oscillator, :Custom,
            :Tensor, :SpinOrbit,
        )
        @test name in names(FewBodyVMC)
        @test getproperty(FewBodyVMC, name) === getproperty(FewBodyHamiltonians, name)
    end
    @test_throws ArgumentError Operators(FewBodyHamiltonians.Hamiltonian())

    masses = [1836.15267343, 3670.48296788, 1.0]
    charges = [1.0, 1.0, -1.0]
    H = Kinetic(masses) + Coulomb(charges)
    legacy = Operators(masses, charges)
    legacy += "Kinetic"
    legacy += "Coulomb"
    lowered = Operators(H)

    @test Matrix(lowered[1].K) == Matrix(legacy[1].K)
    @test coulomb_weights(H) == coulomb_weights(legacy)
    @test [term.coefficient for term in lowered.terms[2:end]] == [1.0, -1.0, -1.0]
    positions = [-1.0 0.0 0.0; 1.0 0.0 0.0; 0.0 0.2 0.0]
    coordinates = jacobi_coordinates(H, positions)
    @test coordinates == jacobi_coordinates(legacy, positions)
    @test particle_coordinates(H, coordinates) == particle_coordinates(legacy, coordinates)
    @test jacobi_coordinates(H, particle_coordinates(H, coordinates)) ≈ coordinates atol = 1.0e-12
    @test_throws DimensionMismatch jacobi_coordinates(H, zeros(2, 3))
    @test_throws DimensionMismatch particle_coordinates(H, zeros(3, 3))

    ψ = DiatomicOneElectronWavefunction(H; α = 1.2, κ = 6.0, R₀ = 2.1)
    legacy_ψ = DiatomicOneElectronWavefunction(legacy; α = 1.2, κ = 6.0, R₀ = 2.1)
    initial = vec(coordinates)
    @test ψ(initial) == legacy_ψ(initial)
    @test local_energy(H, ψ, initial) == local_energy(legacy, legacy_ψ, initial)
    @test local_energy(H, x -> ψ(x), initial) ≈ local_energy(H, ψ, initial) atol = 1.0e-12

    method = VariationalMonteCarlo(
        n_steps = 100, burn_in = 20, n_walkers = 2,
        δ = [12.0, 1.6], r₀ = initial, block_size = 10,
    )
    shared_solution = solve(H, ψ, method)
    legacy_solution = solve(legacy, legacy_ψ, method)
    @test shared_solution.samples == legacy_solution.samples
    @test shared_solution.local_energies == legacy_solution.local_energies
    @test shared_solution.E₀ == legacy_solution.E₀
    @test shared_solution.standard_error == legacy_solution.standard_error
    @test shared_solution.acceptance_rate == legacy_solution.acceptance_rate

    @testset "unsupported systems fail before evaluating the wavefunction" begin
        unevaluated = _ -> error("wavefunction must not be evaluated")
        for invalid in (
                Hamiltonian(Coulomb(1, 2)),
                Kinetic([1.0, 1.0]; units = :nuclear) + Coulomb(1, 2),
                Kinetic([1.0, 1.0]) + Coulomb([1.0, 1.0, -1.0]),
                Kinetic([1.0, 1.0]) + Coulomb(1, 3),
                Kinetic([1.0, 1.0]) + Gaussian(1, 2; exponent = 1.0),
                Kinetic([1.0, 1.0]) + ManyBodyGaussian([1.0 -1.0; -1.0 1.0]),
                Kinetic([1.0, 1.0]) + Oscillator(1, 2),
                Kinetic([1.0, 1.0]) + Custom(r -> r^2, 1, 2),
                Kinetic([1.0, 1.0]) + Tensor(1, 2; exponent = 1.0),
                Kinetic([1.0, 1.0]) + SpinOrbit(1, 2; exponent = 1.0),
            )
            @test_throws ArgumentError Operators(invalid)
            @test_throws ArgumentError solve(invalid, unevaluated, method)
            @test_throws ArgumentError local_energy(invalid, unevaluated, initial)
        end
    end
end
