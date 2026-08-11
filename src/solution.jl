"""Result of a variational Monte Carlo calculation."""
struct VMCSolution{T<:Real,M,L,S}
    E₀::T
    standard_error::T
    variance::T
    acceptance_rate::T
    n_accepted::Int
    n_attempted::Int
    n_samples::Int
    n_discarded::Int
    method::M
    local_energies::L
    samples::S
end

function Base.getproperty(solution::VMCSolution, name::Symbol)
    name === :E && return getfield(solution, :E₀)
    return getfield(solution, name)
end

Base.propertynames(solution::VMCSolution, private::Bool = false) =
    (fieldnames(typeof(solution))..., :E)

"""Return the one-element vector containing the VMC ground-state estimate."""
energies(solution::VMCSolution) = [solution.E₀]

function Base.show(io::IO, solution::VMCSolution)
    print(
        io,
        "VMCSolution(E₀ = $(solution.E₀) ± $(solution.standard_error) Ha, " *
        "acceptance = $(round(solution.acceptance_rate; digits=3)))",
    )
end
