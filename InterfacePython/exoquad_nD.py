# Multi-dimensional experiments for the exotic quadrature paper.
import ExoquadUtils
import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

# ================================================================================================================================
# Experiment #1, f(x) = exp(-|x|^2), rho(x) = sinc(10 * x)
# ================================================================================================================================
pWeightFn = lambda x : np.sinc((10 * x) / np.pi)
pIntegrandFn = lambda xVec : np.exp(-np.linalg.norm(xVec) ** 2)
fBaseIntegral = 0.32099682841103033
for iDim in range(2, 5, 2):
    sTitle = r'Relative Errors for Dimension $n = ' + str(iDim) + r'$'
    sFname = 'exp1_dim' + str(iDim) + '.svg'
    print('Creating ' + sFname + '...')
    ExoquadUtils.plotAccuracy(sTitle, 2000 * (5 ** (iDim - 2)), 1e-14, iDim, fBaseIntegral ** iDim, pIntegrandFn, pWeightFn, 1.0,
                              sFname=sFname)

# ================================================================================================================================
# Experiment #2, f(x) = exp(-Σ xi), rho(x) = sinc(10 * x)
# ================================================================================================================================
pWeightFn = lambda x : np.sinc((10 * x) / np.pi)
pIntegrandFn = lambda xVec : np.exp(np.sum(xVec))
fBaseIntegral = 0.34005259521041936
for iDim in range(2, 5, 2):
    sTitle = r'Relative Errors for Dimension $n = ' + str(iDim) + r'$'
    sFname = 'exp2_dim' + str(iDim) + '.svg'
    print('Creating ' + sFname + '...')
    ExoquadUtils.plotAccuracy(sTitle, 2000 * (5 ** (iDim - 2)), 1e-14, iDim, fBaseIntegral ** iDim, pIntegrandFn, pWeightFn, 1.0,
                              sFname=sFname)

# ================================================================================================================================
# Experiment #3, f(x) = Σ (xi + 1) ^ 1.2, rho(x) = sinc(10 * x)
# ================================================================================================================================
pWeightFn = lambda x : np.sinc((10 * x) / np.pi)
pIntegrandFn = lambda xVec : np.sum(np.power(xVec + 1, 1.2))
fIntegral1D = 0.33379512418920954 # Integral of integrand times the weight function for n=1
fIntegralWeight = 0.3316695188437748 # Integral of weight function
for iDim in range(2, 5, 2):
    sTitle = r'Relative Errors for Dimension $n = ' + str(iDim) + r'$'
    sFname = 'exp3_dim' + str(iDim) + '.svg'
    fIntegral = iDim * (fIntegral1D * fIntegralWeight ** (iDim - 1))
    print('Creating ' + sFname + '...')
    ExoquadUtils.plotAccuracy(sTitle, 2000 * (5 ** (iDim - 2)), 1e-14, iDim, fIntegral, pIntegrandFn, pWeightFn, 1.0,
                              sFname=sFname)
