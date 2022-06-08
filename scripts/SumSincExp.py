# Multi-dimensional experiments for the exotic quadrature paper.
import ExoquadUtils
import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

# f(x) = exp(-Σ xi), rho(x) = sinc(10 * x)
def write_plot(a, iDim, fBaseIntegral, dname):
    pWeightFn = lambda x : np.sinc((10 * x) / np.pi)
    pIntegrandFn = lambda xVec : np.exp(np.sum(xVec))
    sTitle = r'Relative Errors for Dimension $n = ' + str(iDim) + r'$'
    sFname = dname + '/exp2_dim' + str(iDim) + '.svg'
    print('Creating ' + sFname + '...')
    ExoquadUtils.plotAccuracy(sTitle, 2000 * (5 ** (iDim - 2)), 1e-14, iDim, fBaseIntegral ** iDim, pIntegrandFn,
                              pWeightFn, 1.0, sFname=sFname)
