raw"""
    VariationalMonteCarlo(; n_steps=10_000, burn_in=1_000, thinning=1,
                            n_walkers=1, δ=0.5, r₀=[1, 0, 0], block_size=50)

Metropolis VMC options. `n_steps` is the retained count per walker. `δ` may
be a scalar, one width per Jacobi coordinate, or one width per flattened
coordinate. A Jacobi-width vector is repeated for the three Cartesian axes.
The reported standard error is estimated from contiguous block means.
"""
struct VariationalMonteCarlo{T<:AbstractFloat,V<:AbstractVector{T}}
    n_steps::Int
    burn_in::Int
    thinning::Int
    n_walkers::Int
    δ::V
    r₀::V
    block_size::Int
end

function VariationalMonteCarlo(;
    n_steps::Int = 10_000,
    burn_in::Int = 1_000,
    thinning::Int = 1,
    n_walkers::Int = 1,
    δ::Union{Real,AbstractVector{<:Real}} = 0.5,
    r₀::AbstractVector{<:Real} = [1.0, 0.0, 0.0],
    block_size::Int = 50,
)
    n_steps > 0 || throw(ArgumentError("n_steps must be positive"))
    burn_in >= 0 || throw(ArgumentError("burn_in must be nonnegative"))
    thinning > 0 || throw(ArgumentError("thinning must be positive"))
    n_walkers > 0 || throw(ArgumentError("n_walkers must be positive"))
    block_size > 0 || throw(ArgumentError("block_size must be positive"))
    widths = δ isa Real ? [Float64(δ)] : Float64.(δ)
    !isempty(widths) && all(x -> isfinite(x) && x > 0, widths) ||
        throw(ArgumentError("δ must contain positive finite widths"))
    position = Float64.(r₀)
    !isempty(position) && all(isfinite, position) ||
        throw(ArgumentError("r₀ must contain finite coordinates"))
    return VariationalMonteCarlo(
        n_steps,
        burn_in,
        thinning,
        n_walkers,
        widths,
        position,
        block_size,
    )
end

function _proposal_widths(method::VariationalMonteCarlo, dimension::Int)
    length(method.δ) == 1 && return fill(only(method.δ), dimension)
    length(method.δ) == dimension && return method.δ
    dimension % 3 == 0 && length(method.δ) == dimension ÷ 3 && return repeat(method.δ, 3)
    throw(DimensionMismatch("δ must have length 1, $(dimension ÷ 3), or $dimension"))
end

function _kinetic_matrix(operators)
    matrices = [Matrix(term.K) for term in operators if _is_kinetic(term)]
    isempty(matrices) &&
        throw(ArgumentError("the Hamiltonian has no supported kinetic term"))
    K = reduce(+, matrices)
    all(size(matrix) == size(K) for matrix in matrices) ||
        throw(DimensionMismatch("kinetic matrices must have matching dimensions"))
    return K
end

function _potential_energy(operators, position)
    energy = 0.0
    for term in operators
        _is_kinetic(term) && continue
        _is_coulomb(term) ||
            throw(ArgumentError("unsupported VMC operator $(typeof(term))"))
        distance = _pair_distance(position, term.w)
        energy += term.coefficient / distance
    end
    return energy
end

function _weighted_laplacian(hessian, K)
    ncoords = size(K, 1)
    size(hessian) == (3ncoords, 3ncoords) ||
        throw(DimensionMismatch("wavefunction dimension and kinetic matrix do not match"))
    result = zero(eltype(hessian))
    for axis = 0:2, a = 1:ncoords, b = 1:ncoords
        result += K[a, b] * hessian[a+axis*ncoords, b+axis*ncoords]
    end
    return result
end

raw"""
    local_energy(operators, wavefunction, position)

Evaluate ``(H\psi)/\psi``. Arbitrary real callable wavefunctions use a
ForwardDiff Hessian; the bundled diatomic trial function uses an equivalent
analytic expression.
"""
function local_energy(operators, wavefunction, position::AbstractVector{<:Real})
    operators = _vmc_operators(operators)
    value = wavefunction(position)
    iszero(value) && throw(DomainError(value, "wavefunction is zero at position"))
    K = _kinetic_matrix(operators)
    hessian = ForwardDiff.hessian(wavefunction, position)
    return -_weighted_laplacian(hessian, K) / value + _potential_energy(operators, position)
end

function _radial_component_energy(operators, position, wavefunction, electron::Int)
    weights = (wavefunction.w_nn, wavefunction.w_ne1, wavefunction.w_ne2)
    vectors = map(w -> _pair_vector(position, w), weights)
    distances = map(norm, vectors)
    derivatives = (
        2wavefunction.κ * (distances[1] - wavefunction.R₀),
        electron == 1 ? wavefunction.α : zero(wavefunction.α),
        electron == 2 ? wavefunction.α : zero(wavefunction.α),
    )
    second_derivatives = (2wavefunction.κ, 0.0, 0.0)
    K = _kinetic_matrix(operators)
    kinetic = 0.0
    for p = 1:3
        kinetic +=
            dot(weights[p], K, weights[p]) *
            (second_derivatives[p] + 2derivatives[p] / distances[p])
        for q = 1:3
            cosine = dot(vectors[p], vectors[q]) / (distances[p] * distances[q])
            kinetic -=
                derivatives[p] * derivatives[q] * dot(weights[p], K, weights[q]) * cosine
        end
    end
    return kinetic + _potential_energy(operators, position)
end

function local_energy(
    operators,
    wavefunction::DiatomicOneElectronWavefunction,
    position::AbstractVector{<:Real},
)
    operators = _vmc_operators(operators)
    R = _pair_distance(position, wavefunction.w_nn)
    r1 = _pair_distance(position, wavefunction.w_ne1)
    r2 = _pair_distance(position, wavefunction.w_ne2)
    common = -wavefunction.κ * (R - wavefunction.R₀)^2
    log1 = common - wavefunction.α * r1
    log2 = common - wavefunction.α * r2
    shift = max(log1, log2)
    weight1 = exp(log1 - shift)
    weight2 = exp(log2 - shift)
    energy1 = _radial_component_energy(operators, position, wavefunction, 1)
    energy2 = _radial_component_energy(operators, position, wavefunction, 2)
    return (weight1 * energy1 + weight2 * energy2) / (weight1 + weight2)
end

function _log_density(wavefunction, position)
    value = abs2(wavefunction(position))
    value isa Real || throw(ArgumentError("abs2(wavefunction(position)) must be real"))
    return log(value)
end

function _metropolis_samples(wavefunction, method, rng)
    dimension = length(method.r₀)
    widths = _proposal_widths(method, dimension)
    total = method.n_walkers * method.n_steps
    samples = Matrix{Float64}(undef, dimension, total)
    accepted = 0
    attempted = 0
    stored = 0
    transitions = method.burn_in + method.n_steps * method.thinning

    for _ = 1:method.n_walkers
        position = copy(method.r₀)
        log_probability = _log_density(wavefunction, position)
        isfinite(log_probability) ||
            throw(ArgumentError("probability density at r₀ must be positive and finite"))
        for step = 1:transitions
            proposal = position .+ widths .* (rand(rng, dimension) .- 0.5)
            proposed_log_probability = _log_density(wavefunction, proposal)
            attempted += 1
            if isfinite(proposed_log_probability) &&
               log(rand(rng)) < min(0.0, proposed_log_probability - log_probability)
                position = proposal
                log_probability = proposed_log_probability
                accepted += 1
            end
            if step > method.burn_in && (step - method.burn_in) % method.thinning == 0
                stored += 1
                samples[:, stored] = position
            end
        end
    end
    return samples, accepted, attempted
end

function _statistics(values, block_size)
    n = length(values)
    mean_value = sum(values) / n
    variance = n == 1 ? 0.0 : sum(abs2(value - mean_value) for value in values) / (n - 1)
    nblocks = n ÷ block_size
    if nblocks >= 2
        blocks = [
            sum(@view values[((i-1)*block_size+1):(i*block_size)]) / block_size for
            i = 1:nblocks
        ]
        block_mean = sum(blocks) / nblocks
        block_variance = sum(abs2(value - block_mean) for value in blocks) / (nblocks - 1)
        standard_error = sqrt(block_variance / nblocks)
    else
        standard_error = sqrt(variance / n)
    end
    return mean_value, variance, standard_error
end

"""
    solve(operators, wavefunction, method::VariationalMonteCarlo; rng=MersenneTwister(123))

Sample `abs2(wavefunction)` and return a [`VMCSolution`](@ref). The default RNG
makes examples and tests reproducible.
"""
function solve(
    operators,
    wavefunction,
    method::VariationalMonteCarlo;
    rng::Random.AbstractRNG = Random.MersenneTwister(123),
)
    operators = _vmc_operators(operators)
    K = _kinetic_matrix(operators)
    length(method.r₀) == 3size(K, 1) ||
        throw(DimensionMismatch("r₀ must have length $(3size(K, 1))"))
    samples, accepted, attempted = _metropolis_samples(wavefunction, method, rng)
    all_energies =
        [local_energy(operators, wavefunction, sample) for sample in eachcol(samples)]
    finite_energies = filter(isfinite, all_energies)
    isempty(finite_energies) && throw(ArgumentError("all local energies are non-finite"))
    energy, variance, standard_error = _statistics(finite_energies, method.block_size)
    return VMCSolution(
        energy,
        standard_error,
        variance,
        accepted / attempted,
        accepted,
        attempted,
        length(finite_energies),
        length(all_energies) - length(finite_energies),
        method,
        all_energies,
        samples,
    )
end
