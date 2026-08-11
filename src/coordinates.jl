raw"""
    jacobi_transform(masses) -> (J, U)

Return the mass-scaled Jacobi transformation `J` and its Moore–Penrose inverse
`U`. Particle positions `R` (one particle per row) are mapped to internal
coordinates by `X = J * R`. The convention is identical to FewBodyECG.jl.
"""
function jacobi_transform(masses::AbstractVector{<:Real})
    length(masses) >= 2 || throw(ArgumentError("at least two masses are required"))
    all(m -> isfinite(m) && m > 0, masses) ||
        throw(ArgumentError("masses must be positive and finite"))

    mass = Float64.(masses)
    n = length(mass)
    J = zeros(Float64, n - 1, n)
    for k = 1:(n-1)
        partial_mass = sum(@view mass[1:k])
        next_mass = mass[k+1]
        reduced_scale = sqrt(next_mass * partial_mass / (next_mass + partial_mass))
        for j = 1:k
            J[k, j] = reduced_scale * mass[j] / partial_mass
        end
        J[k, k+1] = -reduced_scale
    end
    return J, pinv(J)
end

raw"""
    Λ(masses)

Return the internal kinetic-energy matrix
``\Lambda = J\,\mathrm{diag}(1/(2m_i))\,J^T`` in atomic units.
"""
function Λ(masses::AbstractVector{<:Real})
    J, _ = jacobi_transform(masses)
    return Symmetric(J * Diagonal(0.5 ./ Float64.(masses)) * transpose(J))
end
