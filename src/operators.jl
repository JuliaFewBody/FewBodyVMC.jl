const Operator = FewBodyHamiltonians.Operator

raw"""Kinetic energy ``-\sum_{ab}\Lambda_{ab}\nabla_a\cdot\nabla_b`` in Jacobi coordinates."""
struct KineticOperator{T<:Real,M<:AbstractMatrix{T}} <: FewBodyHamiltonians.KineticTerm
    K::Symmetric{T,M}
end

function KineticOperator(K::AbstractMatrix{<:Real})
    size(K, 1) == size(K, 2) || throw(ArgumentError("K must be square"))
    isapprox(K, transpose(K)) || throw(ArgumentError("K must be symmetric"))
    data = Matrix{Float64}(K)
    isposdef(Symmetric(data)) || throw(ArgumentError("K must be positive definite"))
    return KineticOperator{Float64,typeof(data)}(Symmetric(data))
end

KineticOperator(masses::AbstractVector{<:Real}) = KineticOperator(Λ(masses))

"""Coulomb interaction `coefficient / norm(X' * w)` in Jacobi coordinates."""
struct CoulombOperator{T<:Real,V<:AbstractVector{T}} <: FewBodyHamiltonians.PotentialTerm
    coefficient::T
    w::V
end

function CoulombOperator(coefficient::Real, w::AbstractVector{<:Real})
    isempty(w) && throw(ArgumentError("w must not be empty"))
    all(isfinite, w) || throw(ArgumentError("w must be finite"))
    return CoulombOperator(Float64(coefficient), Float64.(w))
end

"""
    Operators([masses[, charges]])

Hamiltonian builder using the same particle ordering and shorthand as
FewBodyECG.jl. `ops += "Kinetic"` adds the internal kinetic energy and
`ops += "Coulomb"` adds every pair with coefficient `q[i] * q[j]`.
"""
mutable struct Operators
    terms::Vector{Operator}
    masses::Union{Nothing,Vector{Float64}}
    charges::Union{Nothing,Vector{Float64}}
    J::Union{Nothing,Matrix{Float64}}
    U::Union{Nothing,Matrix{Float64}}
end

Operators() = Operators(Operator[], nothing, nothing, nothing, nothing)

function Operators(masses::AbstractVector{<:Real})
    J, U = jacobi_transform(masses)
    return Operators(Operator[], Float64.(masses), nothing, J, U)
end

function Operators(masses::AbstractVector{<:Real}, charges::AbstractVector{<:Real})
    length(masses) == length(charges) ||
        throw(ArgumentError("masses and charges must have the same length"))
    J, U = jacobi_transform(masses)
    return Operators(Operator[], Float64.(masses), Float64.(charges), J, U)
end

function Base.:+(ops::Operators, term::Operator)
    push!(ops.terms, term)
    return ops
end

function Base.:+(ops::Operators, name::AbstractString)
    if name == "Kinetic"
        isnothing(ops.masses) &&
            throw(ArgumentError("\"Kinetic\" requires Operators(masses)"))
        push!(ops.terms, KineticOperator(ops.masses))
    elseif name == "Coulomb"
        isnothing(ops.charges) &&
            throw(ArgumentError("\"Coulomb\" requires Operators(masses, charges)"))
        n = length(ops.masses)
        for i = 1:n, j = (i+1):n
            selector = zeros(Float64, n)
            selector[i] = 1
            selector[j] = -1
            w = transpose(ops.U) * selector
            push!(ops.terms, CoulombOperator(ops.charges[i] * ops.charges[j], w))
        end
    else
        throw(
            ArgumentError(
                "unknown operator \"$name\"; supported: \"Kinetic\", \"Coulomb\"",
            ),
        )
    end
    return ops
end

Base.length(ops::Operators) = length(ops.terms)
Base.getindex(ops::Operators, i::Int) = ops.terms[i]
Base.iterate(ops::Operators, state...) = iterate(ops.terms, state...)

"""Return Coulomb pair weights in insertion order."""
coulomb_weights(ops) = [term.w for term in ops if _is_coulomb(term)]

"""Transform an `N × 3` particle-position matrix to Jacobi coordinates."""
function jacobi_coordinates(ops::Operators, positions::AbstractMatrix{<:Real})
    isnothing(ops.J) && throw(ArgumentError("Operators does not contain masses"))
    size(positions) == (length(ops.masses), 3) ||
        throw(DimensionMismatch("positions must have size ($(length(ops.masses)), 3)"))
    return ops.J * positions
end

"""Back-transform Jacobi coordinates to centre-of-mass-free particle positions."""
function particle_coordinates(ops::Operators, coordinates::AbstractMatrix{<:Real})
    isnothing(ops.U) && throw(ArgumentError("Operators does not contain masses"))
    size(coordinates) == (length(ops.masses) - 1, 3) || throw(
        DimensionMismatch("coordinates must have size ($(length(ops.masses) - 1), 3)"),
    )
    return ops.U * coordinates
end

_is_kinetic(term) = term isa FewBodyHamiltonians.KineticTerm && hasproperty(term, :K)
_is_coulomb(term) =
    term isa FewBodyHamiltonians.PotentialTerm &&
    nameof(typeof(term)) == :CoulombOperator &&
    hasproperty(term, :coefficient) &&
    hasproperty(term, :w)

function Base.show(io::IO, ops::Operators)
    suffix = length(ops) == 1 ? "" : "s"
    print(io, "Operators with $(length(ops)) term$suffix")
end
