# Basic plots for the exotic quadrature paper.
import TasmanianSG
import ExoquadUtils
import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

# ================================================================================================================================
# Basic plots for the integral of sinc(a*x) over [-1,1].
# ================================================================================================================================

# Parameters
iMaxPoints = 30
fMinErr = 1e-13

# Initialization.
iFontSize = 15
plt.rc('axes', titlesize=iFontSize)
plt.rc('axes', labelsize=iFontSize)
plt.rc('xtick', labelsize=iFontSize-2)
plt.rc('ytick', labelsize=iFontSize-2)
plt.rc('figure', titlesize=iFontSize+2)
plt.figure()

for (a, m) in [(1,'o'), (5,'v'), (10,'v'), (25,'D')]:
    fIntegral = 2 * sp.sici(a)[0] / a
    pIntegrandFn = lambda x : np.sinc(a * x / np.pi)
    print('Creating plots for a=' + str(a))
    fErrLower = 1e-20
    pGridFn = lambda iDepth : TasmanianSG.makeGlobalGrid(1, 0, iDepth, 'qptotal', 'gauss-legendre')
    lNumPoints, lErr, _ = ExoquadUtils.getData(pGridFn, iMaxPoints, fMinErr, fIntegral, pIntegrandFn)
    plt.plot(lNumPoints, fErrLower + lErr, label=r'$a=' + str(a) + r'$', marker=m)

# Show plots.
plt.title("Relative Errors for Different Frequencies")
plt.xlabel("Number of Points")
plt.ylabel("Relative Error")
plt.yscale("log")
plt.grid(linestyle="--", linewidth=0.25)
plt.legend()
plt.savefig('sinc_plot.svg')

