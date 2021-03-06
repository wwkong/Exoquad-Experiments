# Multi-dimensional experiments for the exotic quadrature paper.
import ExoquadUtils
import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

# f(x) = Σ (xi + 1) ^ 1.2, rho(x) = sinc(20 * x)
# fIntegral1D = Integral of integrand times the weight function for n=1
# fIntegralWeight = Integral of weight function
def write_plot(a, iDim, fIntegral1D, fIntegralWeight, dname):
    sFname = dname + '/exp3_dim' + str(iDim) + '.svg'
    print('Creating ' + sFname + '...', end='', flush=True)
    pWeightFn = lambda x : np.sinc((a * x) / np.pi)
    pIntegrandFn = lambda xVec : np.sum(np.power(xVec + 1, 1.2))
    sTitle = r'Relative Errors for Dimension $n = ' + str(iDim) + r'$'
    fIntegral = iDim * (fIntegral1D * fIntegralWeight ** (iDim - 1))
    ExoquadUtils.plotAccuracy(sTitle, 2000 * (5 ** (iDim - 2)), 1e-10, iDim, fIntegral, pIntegrandFn,
                              pWeightFn, 1.0, sFname=sFname, iNref=500*iDim)
    print("DONE")
