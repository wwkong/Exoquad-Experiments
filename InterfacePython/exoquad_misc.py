# Basic plots for the exotic quadrature paper.
import TasmanianSG
import ExoquadUtils
import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

# ================================================================================================================================
# Experiment #1, f(x) = exp(-||x||^2), rho(x) = sinc(20 * x)
# ================================================================================================================================
pWeightFn = lambda x : np.sinc((20 * x) / np.pi)
pIntegrandFn = lambda xVec : np.exp(-np.linalg.norm(xVec) ** 2)
fBaseIntegral = 0.15609592912582969
for iDim in [1, 2]:
    sTitle = r'Relative Errors for Dimension $n = ' + str(iDim) + r'$'
    sFname = 'exp1_dim' + str(iDim) + '.svg'
    print('Creating ' + sFname + '...')
    ExoquadUtils.plotAccuracy(sTitle, 2000 * (5 ** (iDim - 2)), 1e-12, iDim, fBaseIntegral ** iDim, pIntegrandFn, pWeightFn, 1.0,
                              sFname=sFname)

# # ================================================================================================================================
# # Basic plots for the integral of sinc(a*x) over [-1,1].
# # ================================================================================================================================

# # Parameters
# iMaxPoints = 30
# fMinErr = 1e-13

# # Initialization.
# iFontSize = 15
# plt.rc('axes', titlesize=iFontSize)
# plt.rc('axes', labelsize=iFontSize)
# plt.rc('xtick', labelsize=iFontSize-2)
# plt.rc('ytick', labelsize=iFontSize-2)
# plt.rc('figure', titlesize=iFontSize+2)
# plt.figure()

# for (a, m) in [(1,'o'), (5,'v'), (10,'v'), (25,'D')]:
#     fIntegral = 2 * sp.sici(a)[0] / a
#     pIntegrandFn = lambda x : np.sinc(a * x / np.pi)
#     print('Creating plots for a=' + str(a))
#     fErrLower = 1e-20
#     pGridFn = lambda iDepth : TasmanianSG.makeGlobalGrid(1, 0, iDepth, 'qptotal', 'gauss-legendre')
#     lNumPoints, lErr, _ = ExoquadUtils.getData(pGridFn, iMaxPoints, fMinErr, fIntegral, pIntegrandFn)
#     plt.plot(lNumPoints, fErrLower + lErr, label=r'$a=' + str(a) + r'$', marker=m)

# # Show plots.
# plt.title("Relative Errors for Different Frequencies")
# plt.xlabel("Number of Points")
# plt.ylabel("Relative Error")
# plt.yscale("log")
# plt.grid(linestyle="--", linewidth=0.25)
# plt.legend()
# plt.savefig('sinc_plot.svg')

