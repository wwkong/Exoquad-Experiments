using GaussQuadrature
using Random
using SpecialFunctions

# Create a quadrature for ∫ μ(τ) * exp(im * c * τ * cos(θ)) dτ where:
#
#     μ(x) = I₀(π * sqrt(1 - x²)) where I₀ is the modified Bessel function of order zero.
#

# Initialize.
a = 10 * pi
N = 300
aux_int = y -> sinc(sqrt(Complex(y ^ 2 - pi ^ 2)) / pi)
f = x -> exp(im * a * x)
sample = a .* range(-N, N, step=1) ./ N
moments = aux_int.(sample)
integralVal = aux_int.(a)
#display(integralVal)

# Evaluate.
num_points = []
relative_err = []
tols = [10 .^ (-0.5:-0.1:-0.9) ; 10.0 .^ (-1:-1:-14)]
for tol in tols
    nodes, weights = GaussQuadrature.ExpGauss.getExpGaussNodesAndWeights2(a, moments, tol)
    fvals = f.(nodes)
    integralApprox = sum(fvals .* weights) 
    append!(num_points, length(weights))
    append!(relative_err, abs(integralApprox - integralVal) / abs(integralVal))
end

display(hcat(num_points, relative_err))
