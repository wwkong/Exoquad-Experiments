__precompile__()
module ExpGauss
#=
Implementation of the algorithms in the papers:

[1] Beylkin, G., & Monzon, L. (2002). On generalized Gaussian quadratures for exponentials and their applications. Applied and
Computational Harmonic Analysis, 12(3), 332-373.

[2] Reynolds, M., Beylkin, G., & Monzón, L. (2013). On generalized Gaussian quadratures for bandlimited exponentials. Applied and
Computational Harmonic Analysis, 34(3), 352-365.

=#

using JuMP
using ECOS
using ToeplitzMatrices: Toeplitz
using LinearAlgebra: eigen, svd, pinv
using Polynomials: Polynomial, roots
using SpecialMatrices: Vandermonde

function getExpGaussNodesAndWeights1(c, moments, tol)
    # Based on the quadrature algorithm in [1, Algorithm 2]. Computes quadrature for integrals of the form
    #
    #   I(x,c) = ∫ exp(i*x*c*τ) μ(τ) dτ   for every  x ∈ [-1,1],
    #
    # using the exponential moments of μ. In particular, if `moments` is length 2*N+1, then this method expects
    #
    #   moments[i] = ∫ exp(i*c*τ*[i-N-1]/N) μ(τ) dτ
    #
    if length(moments) % 2 != 1 || length(moments) <= 3
        throw(ArgumentError("Value array must be an odd length that is greater than 3!"))
    end
    N = Int((length(moments) - 1) / 2)
    # Generate the nodes.
    c1 = moments[(N+1):-1:1]
    r1 = moments[(N+1):1:end]
    EigTn = eigen(Matrix(Toeplitz(c1, r1)))
    idx = argmin(abs.(EigTn.values .- tol))
    coeffs = EigTn.vectors[:, idx]
    gamma_unsrt = roots(Polynomial(coeffs))
    radii = abs.(gamma_unsrt)
    nodes = sort!(imag.(log.(gamma_unsrt ./ radii)) ./ c .* N)
    gamma = exp.(im .* nodes .* c ./ N) .* radii
    # Generate the weights based on the Vandermonde system in [1, Eq. (4.18)].
    M = length(gamma)
    V = zeros(ComplexF64, (M, M))
    for i=1:M
        for j=1:M
            V[i, j] = gamma[j] ^ (-i)
        end
    end
    b = moments[N:-1:(N-M+1)]
    weights = V \ b
    # Output.
    return nodes, weights
end

function getExpGaussNodesAndWeights2(c, moments, tol)
    # Based on the quadrature algorithm in [2, Algorithm 1]. Expects the same inputs as in getExpGaussNodesAndWeights1().
    if length(moments) % 2 != 1 || length(moments) <= 3
        throw(ArgumentError("Value array must be an odd length that is greater than 3!"))
    end
    N = Int((length(moments) - 1) / 2)
    # Generate the nodes.
    c1 = moments[(N+1):-1:1]
    r1 = moments[(N+1):1:end]
    SvdTn = svd(Matrix(Toeplitz(c1, r1)))
    idx = argmin(abs.(SvdTn.S ./ SvdTn.S[1] .- tol))
    UM_tilde = SvdTn.U[1:N, 1:idx-1]
    UM_hat = SvdTn.U[2:N+1, 1:idx-1]
    CM = pinv(UM_tilde) * UM_hat;
    gamma_unsrt = eigen(CM).values
    radii = abs.(gamma_unsrt)
    nodes = sort!(imag.(log.(gamma_unsrt ./ radii)) ./ c .* N)
    gamma = exp.(im .* nodes .* c ./ N) .* radii
    # Generate the weights by ℓ∞-minimization of the system [V * `weights` = `moments`].
    M = length(gamma)
    V = zeros(ComplexF64, (2*N+1, M))
    for i=-N:1:N
        for j=1:M
            V[i+N+1, j] = gamma[j] ^ (i)
        end
    end
    # Generate the minimizer by solver an SOCP described in [2, Subsection 8.2.1].
    SOCP = JuMP.Model(ECOS.Optimizer)
    JuMP.set_silent(SOCP)
    @variable(SOCP, re_weights[1:M])
    @variable(SOCP, im_weights[1:M])
    @variable(SOCP, t)
    for i=1:2*N+1
        a = V[i,:]'
        A = vcat(hcat(real.(a), -imag.(a)),
                  hcat(imag.(a),  real.(a)))
        b = [real(moments[i]); imag(moments[i])]
        @constraint(SOCP, vcat(t, A * [re_weights; im_weights] .- b) in SecondOrderCone())
    end
    @objective(SOCP, Min, t)
    JuMP.optimize!(SOCP)
    if (JuMP.termination_status(SOCP) == MOI.OPTIMAL ||
        JuMP.termination_status(SOCP) == MOI.ALMOST_OPTIMAL ||
        JuMP.termination_status(SOCP) == MOI.TIME_LIMIT && JuMP.has_values(SOCP))
        weights = JuMP.value.(re_weights) .+ im .* JuMP.value.(im_weights)
    else
        error("The model terminated improperly with status: " * string(JuMP.termination_status(SOCP)))
    end
    # Output.
    return nodes, weights
end

end # module
