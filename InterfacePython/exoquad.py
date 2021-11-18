# Experiments for the exotic quadrature paper.
import TasmanianSG
import TasmanianAddons

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['text.usetex'] = True

def getData(pGridConstructor, iMaxPoints, fMinErr, fIntegral, pIntegrandFn):
# Given a grid constructor with single input iDepth, returns the # of points and accuracy vectors of the grids up to a given
# maximum number of points.
    iDepth = 0
    iNumPoints = 0
    lNumPoints = np.array([])
    fErr = float("inf")
    lErr = np.array([])
    while (iNumPoints < iMaxPoints and fErr > fMinErr):
        pGrid = pGridConstructor(iDepth)
        points = pGrid.getPoints()
        weights = pGrid.getQuadratureWeights()
        fVals = np.apply_along_axis(pIntegrandFn, 1, points)
        fErr = abs(fIntegral - np.dot(fVals, weights)) / abs(fIntegral)
        lErr = np.append(lErr, fErr)
        iNumPoints = pGrid.getNumPoints()
        lNumPoints = np.append(lNumPoints, iNumPoints)
        iDepth = iDepth + 1
    return lNumPoints, lErr

def plotAccuracy(sTitle, iMaxPoints, fMinErr, iDim, fIntegral, pIntegrandFn, pWeightFn, fShift, sFname = "plot.svg", iSymmetric = False, iNref = 200):
# For an upper bound on the number of points and a given dimension, plot the accuracy of Gauss-Legendre, Gauss-Legendre (odd),
# Gauss-Patterson, and Exotic Quadrature grids with respect to a given integral value fIntegral, a 1D weight function pWeightFn,
# and an nD integrand function pIntegrandFn.

    # Initialization.
    iFontSize = 15
    plt.rc('axes', titlesize=iFontSize)
    plt.rc('axes', labelsize=iFontSize)
    plt.rc('xtick', labelsize=iFontSize-2)
    plt.rc('ytick', labelsize=iFontSize-2)
    plt.rc('figure', titlesize=iFontSize+2)
    plt.figure()
    pCombinedFn = lambda xVec : pIntegrandFn(xVec) * np.prod(pWeightFn(xVec))
    fErrLower = 1e-20

    # Gauss-Legendre (all).
    pGridFn = lambda iDepth : TasmanianSG.makeGlobalGrid(iDim, 0, iDepth, 'qptotal', 'gauss-legendre')
    lNumPoints, lErr = getData(pGridFn, iMaxPoints, fMinErr, fIntegral, pCombinedFn)
    plt.plot(lNumPoints,  + lErr, label="GL", marker='o')

    # Gauss-Legendre (odd).
    pGridFn = lambda iDepth : TasmanianSG.makeGlobalGrid(iDim, 0, iDepth, 'qptotal', 'gauss-legendre-odd')
    lNumPoints, lErr = getData(pGridFn, iMaxPoints, fMinErr, fIntegral, pCombinedFn)
    plt.plot(lNumPoints, fErrLower + lErr, label="GL-odd", marker='v')

    # Gauss-Patterson.
    pGridFn = lambda iDepth : TasmanianSG.makeGlobalGrid(iDim, 0, iDepth, 'qptotal', 'gauss-patterson')
    lNumPoints, lErr = getData(pGridFn, iMaxPoints, fMinErr, fIntegral, pCombinedFn)
    plt.plot(lNumPoints, fErrLower + lErr, label="GP", marker='s')

    # Gauss-Exotic.
    def pExoGridFn(iDepth):
        exoCT = TasmanianAddons.createExoticQuadratureFromFunction(iDepth + 1, fShift, pWeightFn, iNref, "Exotic Quadrature", iSymmetric)
        return TasmanianSG.makeGlobalGridCustom(iDim, 0, iDepth, 'qptotal', exoCT)
    lNumPoints, lErr = getData(pExoGridFn, iMaxPoints, fMinErr, fIntegral, pIntegrandFn)
    plt.plot(lNumPoints, fErrLower + lErr, label="Exotic", marker='*')

    # Show plots.
    plt.title(sTitle)
    plt.xlabel("Number of Points")
    plt.ylabel("Relative Error")
    plt.yscale("log")
    plt.grid(linestyle="--", linewidth=0.25)
    plt.legend(loc="upper right")
    plt.savefig(sFname)

# ================================================================================================================================
# # Unit Test
# ================================================================================================================================
# sincFn = lambda x : np.sinc(10 * x / np.pi)
# expNegSqrNormFn = lambda xVec : np.exp(-np.linalg.norm(xVec) ** 2)
# iMaxPoints = 5000
# iDim = 4
# fIntegral = 0.32099682841103033 ** iDim
# fShift = 1.0
# plotAccuracy("Test", iMaxPoints, iDim, fIntegral, expNegSqrNormFn, sincFn, fShift)

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
    plotAccuracy(sTitle, 2000 * (5 ** (iDim - 2)), 1e-14, iDim, fBaseIntegral ** iDim, pIntegrandFn, pWeightFn, 1.0, sFname=sFname)

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
    plotAccuracy(sTitle, 2000 * (5 ** (iDim - 2)), 1e-14, iDim, fBaseIntegral ** iDim, pIntegrandFn, pWeightFn, 1.0, sFname=sFname)

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
    plotAccuracy(sTitle, 2000 * (5 ** (iDim - 2)), 1e-14, iDim, fIntegral, pIntegrandFn, pWeightFn, 1.0, sFname=sFname)


