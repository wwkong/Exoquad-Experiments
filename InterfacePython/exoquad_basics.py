# Basic plots for the exotic quadrature paper.
import TasmanianSG
import TasmanianAddons

import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp

plt.rcParams['text.usetex'] = True

def getData(pGridConstructor, iMaxPoints, fMinErr, fIntegral, pIntegrandFn):
# Given a grid constructor with single input iDepth, returns the # of points and accuracy vectors of the grids up to a given
# maximum number of points.
    iDepth = 0
    lEstimates = np.array([])
    iNumPoints = 0
    lNumPoints = np.array([])
    fErr = float("inf")
    lErr = np.array([])
    while (iNumPoints < iMaxPoints and fErr > fMinErr):
        pGrid = pGridConstructor(iDepth)
        points = pGrid.getPoints()
        weights = pGrid.getQuadratureWeights()
        fVals = np.squeeze(np.apply_along_axis(pIntegrandFn, 1, points))
        fEstimate = np.dot(fVals, weights)
        lEstimates = np.append(lEstimates, fEstimate)
        fErr = abs(fIntegral - fEstimate) / abs(fIntegral)
        lErr = np.append(lErr, fErr)
        iNumPoints = pGrid.getNumPoints()
        lNumPoints = np.append(lNumPoints, iNumPoints)
        iDepth = iDepth + 1
    return lNumPoints, lErr, lEstimates

# ================================================================================================================================
# Basic plots for the integral of sinc(a*x) over [-1,1].
# ================================================================================================================================

# Parameters
iMaxPoints = 30
fMinErr = -1

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
    print(fIntegral)
    fErrLower = 1e-20
    pGridFn = lambda iDepth : TasmanianSG.makeGlobalGrid(1, 0, iDepth, 'qptotal', 'gauss-legendre')
    lNumPoints, lErr, _ = getData(pGridFn, iMaxPoints, fMinErr, fIntegral, pIntegrandFn)
    plt.plot(lNumPoints, fErrLower + lErr, label=r'$a=' + str(a) + r'$', marker=m)

# Show plots.
plt.title("Relative Errors for Different Frequencies")
plt.xlabel("Number of Points")
plt.ylabel("Relative Error")
plt.yscale("log")
plt.grid(linestyle="--", linewidth=0.25)
plt.legend()
plt.savefig('sinc_plot.svg')
