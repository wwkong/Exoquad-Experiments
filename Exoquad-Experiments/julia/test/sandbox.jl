using GaussQuadrature
using Random

# --------------------------------------------------------------------------------------------------------------------------------
# EXAMPLE 1.
# --------------------------------------------------------------------------------------------------------------------------------

# Create a quadrature for ∫ μ(x) * exp(im * c * τ) dτ where:
#
#     μ(x) = 1.0 ? (-1/6 <= x && x <= 1/6) : 0.0
#
c = 15 * π                                    # Bandwidth limit.
N = 97                                        # Sampling frequency; should be higher than the Nyquist rate of 2c/π.
samples = range(-N, N, step=1) ./ N           # Sample points k/N of μ(x) for k=-N,...,N.
moments = sinc.(c .* samples ./ (6 * π)) ./ 3 # Values of ∫ μ(τ) * exp(im * c * τ * k/N) dτ for k=-N,...,N.

# # Generate the nodes and weights.
# tol = 1e-12
# nodes, weights = GaussQuadrature.ExpGauss.getExpGaussNodesAndWeights1(c, moments, tol)
# display(hcat(nodes, weights))

# --------------------------------------------------------------------------------------------------------------------------------
# EXAMPLE 2.
# --------------------------------------------------------------------------------------------------------------------------------
# Create a quadrature for ∫ μ(x) * exp(im * c * τ) dτ where:
#
#     μ(x) = I₀(π * sqrt(1-x²)) where I₀ is the modified Bessel function of order zero.
#
c = 10 * pi
N = 252
sample = c .* range(-N, N, step=1) ./ N
moments = sinc.(sqrt.(Complex.(sample.^2 .- pi^2)) ./ pi)

# Generate the nodes and weights.
tol = 1.2 * 10^(-15)
nodes, weights = GaussQuadrature.ExpGauss.getExpGaussNodesAndWeights1(c, moments, tol)
display(hcat(nodes, weights))
println("\n")
nodes, weights = GaussQuadrature.ExpGauss.getExpGaussNodesAndWeights2(c, moments, tol)
display(hcat(nodes, weights))
