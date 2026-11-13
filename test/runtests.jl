using FewBodyVMC
using LinearAlgebra
using Test

include("shared_hamiltonians.jl")

function hdplus_problem()
    proton_mass = 1836.15267343
    deuteron_mass = 3670.48296788
    ops = Operators([proton_mass, deuteron_mass, 1.0], [1.0, 1.0, -1.0])
    ops += "Kinetic"
    ops += "Coulomb"
    wavefunction = DiatomicOneElectronWavefunction(ops; α = 1.2, κ = 6.0, R₀ = 2.1)
    positions = [-1.0 0.0 0.0; 1.0 0.0 0.0; 0.0 0.2 0.0]
    initial = vec(jacobi_coordinates(ops, positions))
    return ops, wavefunction, initial
end

@testset "FewBodyVMC.jl" begin
    @testset "Jacobi coordinates and operators" begin
        masses = [1836.15267343, 3670.48296788, 1.0]
        charges = [1.0, 1.0, -1.0]
        J, U = jacobi_transform(masses)

        @test size(J) == (2, 3)
        @test size(U) == (3, 2)
        @test J * U ≈ I(2) atol = 1.0e-12
        @test Matrix(Λ(masses)) ≈ 0.5I(2) atol = 1.0e-12

        ops = Operators(masses, charges)
        ops += "Kinetic"
        ops += "Coulomb"
        @test length(ops) == 4
        @test length(coulomb_weights(ops)) == 3
        @test_throws ArgumentError ops + "unsupported"

        positions = [-1.0 0.0 0.0; 1.0 0.0 0.0; 0.0 0.2 0.0]
        coordinates = jacobi_coordinates(ops, positions)
        @test jacobi_coordinates(ops, particle_coordinates(ops, coordinates)) ≈ coordinates atol =
            1.0e-12
    end

    @testset "local energy" begin
        kinetic = Operators()
        kinetic += KineticOperator(reshape([0.5], 1, 1))
        ψ(x) = exp(-sum(abs2, x) / 2)
        position = [0.2, -0.3, 0.4]
        @test local_energy(kinetic, ψ, position) ≈ 1.5 - sum(abs2, position) / 2 atol =
            1.0e-12

        ops, wavefunction, initial = hdplus_problem()
        analytic = local_energy(ops, wavefunction, initial)
        autodiff = local_energy(ops, x -> wavefunction(x), initial)
        @test analytic ≈ autodiff atol = 1.0e-12
    end

    @testset "HD+ ground state" begin
        ops, wavefunction, initial = hdplus_problem()
        method = VariationalMonteCarlo(
            n_steps = 2_000,
            burn_in = 500,
            thinning = 2,
            n_walkers = 2,
            δ = [12.0, 1.6],
            r₀ = initial,
            block_size = 50,
        )
        solution = solve(ops, wavefunction, method)

        @test -0.65 < solution.E₀ < -0.5
        @test solution.E == solution.E₀
        @test energies(solution) == [solution.E₀]
        @test solution.n_samples == method.n_steps * method.n_walkers
        @test solution.n_discarded == 0
        @test size(solution.samples) == (6, solution.n_samples)
        @test 0.3 < solution.acceptance_rate < 0.8
        @test isfinite(solution.standard_error)
    end

    @testset "validation" begin
        @test_throws ArgumentError VariationalMonteCarlo(n_steps = 0)
        @test_throws ArgumentError VariationalMonteCarlo(burn_in = -1)
        @test_throws ArgumentError VariationalMonteCarlo(δ = 0)
        @test_throws ArgumentError jacobi_transform([1.0])
        @test_throws ArgumentError jacobi_transform([1.0, -1.0])
    end
end
