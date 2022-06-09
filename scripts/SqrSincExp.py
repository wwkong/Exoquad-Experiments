# Multi-dimensional experiments for the exotic quadrature paper.
import ExoquadUtils
import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

# f(x) = exp(-|x|^2), rho(x) = sinc(10 * x)
def write_plot(a, iDim, fBaseIntegral, dname):
    sFname = dname + '/exp1_dim' + str(iDim) + '.svg'
    print('Creating ' + sFname + '...', end='', flush=True)
    pWeightFn = lambda x : np.sinc((a * x) / np.pi)
    pIntegrandFn = lambda xVec : np.exp(-np.linalg.norm(xVec) ** 2)
    sTitle = r'Relative Errors for Dimension $n = ' + str(iDim) + r'$'
    ExoquadUtils.plotAccuracy(sTitle, 2000 * (5 ** (iDim - 2)), 1E-14, iDim, fBaseIntegral ** iDim, pIntegrandFn,
                              pWeightFn, 1.0, sFname=sFname, iNref=500*iDim)
    print("DONE")
