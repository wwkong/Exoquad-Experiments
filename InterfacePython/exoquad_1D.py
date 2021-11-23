# One-dimensional experiments for the exotic quadrature paper.
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

def plotAccuracy(sTitle, iMaxPoints, fMinErr, iDim, fIntegral, pIntegrandFn, pWeightFnReal, pWeightFnIm, fShift,
                 lGBNumPoints, lGBRelativeErr, sFname = "plot.svg", iSymmetric = False, iNref = 300):
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
    pCombinedFn = lambda x : (pWeightFnReal(x) +  pWeightFnIm(x) * 1j) * pIntegrandFn(x)
    fErrLower = 1e-20

    # Gauss-Legendre (all).
    pGridFn = lambda iDepth : TasmanianSG.makeGlobalGrid(iDim, 0, iDepth, 'qptotal', 'gauss-legendre')
    lNumPoints, lErr, _ = getData(pGridFn, iMaxPoints, fMinErr, fIntegral, pCombinedFn)
    plt.plot(lNumPoints, fErrLower + lErr, label="GL", marker='o')

    # Gauss-Legendre (odd).
    pGridFn = lambda iDepth : TasmanianSG.makeGlobalGrid(iDim, 0, iDepth, 'qptotal', 'gauss-legendre-odd')
    lNumPoints, lErr, _ = getData(pGridFn, iMaxPoints, fMinErr, fIntegral, pCombinedFn)
    plt.plot(lNumPoints, fErrLower + lErr, label="GL-odd", marker='v')

    # Gauss-Patterson.
    pGridFn = lambda iDepth : TasmanianSG.makeGlobalGrid(iDim, 0, iDepth, 'qptotal', 'gauss-patterson')
    lNumPoints, lErr, _ = getData(pGridFn, iMaxPoints, fMinErr, fIntegral, pCombinedFn)
    plt.plot(lNumPoints, fErrLower + lErr, label="GP", marker='s')

    # Gauss-Exotic.
    dummy = 1.0
    def pExoGridFn(iDepth, pWeightFn):
        exoCT = TasmanianAddons.createExoticQuadratureFromFunction(iDepth + 1, fShift, pWeightFn, iNref, "Exotic Quadrature", iSymmetric)
        return TasmanianSG.makeGlobalGridCustom(iDim, 0, iDepth, 'qptotal', exoCT)
    pExoGridFnReal = lambda iDepth : pExoGridFn(iDepth, pWeightFnReal)
    lNumPoints, _, lEstReal = getData(pExoGridFnReal, iMaxPoints, fMinErr, dummy, pIntegrandFn)
    pExoGridFnIm = lambda iDepth : pExoGridFn(iDepth, pWeightFnIm)
    _, _, lEstIm = getData(pExoGridFnIm, iMaxPoints, fMinErr, dummy, pIntegrandFn)
    lErr = abs(lEstReal + lEstIm * 1j - fIntegral) / abs(fIntegral)
    plt.plot(lNumPoints, fErrLower + lErr, label="Exotic", marker='*')

    # Gauss-Bandlimited
    plt.plot(lGBNumPoints, fErrLower + lGBRelativeErr, label="GB", marker='D')

    # Show plots.
    plt.title(sTitle)
    plt.xlabel("Number of Points")
    plt.ylabel("Relative Error")
    plt.yscale("log")
    plt.grid(linestyle="--", linewidth=0.25)
    plt.legend(loc="upper right")
    plt.savefig(sFname)

# ================================================================================================================================
# Experiment #1, f(x) = exp(im * a * x), rho(x) = I0(pi * sqrt(1 - x^2)) / 2
# ================================================================================================================================
# Raw arrays from Julia (a=10).
lGBPoints = np.array([8,9,10,12,13,14,15,16,17,18,19,20,21,22])
lGBErrs = np.array([41.0496,56.369,57.0769,3.87772,0.587056,0.0466392,0.023792,0.0265805,0.02104,0.0176328,0.0130911,0.00341601,
                    0.00912594,0.0602533])
# Python arrays (a=10).
a = 10 * np.pi
pIntegrandFn = lambda x : sp.iv(0, np.pi * np.sqrt(1 - x ** 2)) / 2
pWeightFnReal = lambda x : np.cos(a * x)
pWeightFnIm = lambda x : np.sin(a * x)
fBaseIntegral = np.sinc(np.sqrt(a ** 2 - np.pi ** 2) / np.pi)
sTitle = r'Relative Errors for Frequency $a=' + str(int(a / np.pi)) +  r'$'
sFname = 'radio1_a10.svg'
print('Creating ' + sFname + '...')
print(fBaseIntegral)
plotAccuracy(sTitle, 25, 1e-14, 1, fBaseIntegral, pIntegrandFn, pWeightFnReal, pWeightFnIm, 1.0, lGBPoints, lGBErrs, sFname=sFname)

# Raw arrays from Julia (a=20).
lGBPoints = np.array([17,18,19,20,22,23,24,25,27,27,28,29,30,31,32,33,34])
lGBErrs = np.array([85.7029,75.4862,1.32635,137.758,16.786,3.58022,0.441072,0.148161,0.181965,0.181965,0.172947,0.129785,
                    0.111163,0.105705,0.0516552,0.0990842,0.796096])
# Python arrays (a=20).
a = 20 * np.pi
pIntegrandFn = lambda x : sp.iv(0, np.pi * np.sqrt(1 - x ** 2)) / 2
pWeightFnReal = lambda x : np.cos(a * x)
pWeightFnIm = lambda x : np.sin(a * x)
fBaseIntegral = np.sinc(np.sqrt(a ** 2 - np.pi ** 2) / np.pi)
sTitle = r'Relative Errors for Frequency $a=' + str(int(a / np.pi)) + r'$'
sFname = 'radio1_a20.svg'
print('Creating ' + sFname + '...')
print(fBaseIntegral)
plotAccuracy(sTitle, 31, 1e-14, 1, fBaseIntegral, pIntegrandFn, pWeightFnReal, pWeightFnIm, 1.0, lGBPoints, lGBErrs, sFname=sFname)

