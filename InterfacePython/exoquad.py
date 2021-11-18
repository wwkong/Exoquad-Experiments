# Experiments for the exotic quadrature paper.
import TasmanianSG
import TasmanianAddons

import matplotlib.pyplot as plt
import numpy as np

def getData(pGridConstructor, iMaxPoints, fIntegral, pIntegrandFn):
# Given a grid constructor with single input iDepth, returns the # of points and accuracy vectors of the grids up to a given
# maximum number of points.
    iDepth = 0
    iNumPoints = 0
    lNumPoints = np.array([])
    lAcc = np.array([])
    while (iNumPoints < iMaxPoints):
        pGrid = pGridConstructor(iDepth)
        points = pGrid.getPoints()
        weights = pGrid.getQuadratureWeights()
        fVals = np.apply_along_axis(pIntegrandFn, 1, points)
        lAcc = np.append(lAcc, abs(fIntegral - np.dot(fVals, weights)))
        iNumPoints = pGrid.getNumPoints()
        lNumPoints = np.append(lNumPoints, iNumPoints)
        iDepth = iDepth + 1
    return lNumPoints, lAcc

def plotAccuracy(sTitle, iMaxPoints, iDim, fIntegral, pIntegrandFn, pWeightFn, fShift, iSymmetric = False, iNref = 100):
# For an upper bound on the number of points and a given dimension, plot the accuracy of Gauss-Legendre, Gauss-Legendre (odd),
# Gauss-Patterson, and Exotic Quadrature grids with respect to a given integral value fIntegral, a 1D weight function pWeightFn,
# and an nD integrand function pIntegrandFn.

    # Helper variables.
    pCombinedFn = lambda xVec : pIntegrandFn(xVec) * np.prod(pWeightFn(xVec))
    
    # Gauss-Legendre (all).
    pGridFn = lambda iDepth : TasmanianSG.makeGlobalGrid(iDim, 0, iDepth, 'qptotal', 'gauss-legendre')
    lNumPoints, lAcc = getData(pGridFn, iMaxPoints, fIntegral, pCombinedFn)
    plt.plot(lNumPoints, np.log(1e-15 + lAcc), label="Gauss-Legendre (all)", marker='o')

    # Gauss-Legendre (odd).
    pGridFn = lambda iDepth : TasmanianSG.makeGlobalGrid(iDim, 0, iDepth, 'qptotal', 'gauss-legendre-odd')
    lNumPoints, lAcc = getData(pGridFn, iMaxPoints, fIntegral, pCombinedFn)
    plt.plot(lNumPoints, np.log(1e-15 + lAcc), label="Gauss-Legendre (odd)", marker='v')

    # Gauss-Patterson.
    pGridFn = lambda iDepth : TasmanianSG.makeGlobalGrid(iDim, 0, iDepth, 'qptotal', 'gauss-patterson')
    lNumPoints, lAcc = getData(pGridFn, iMaxPoints, fIntegral, pCombinedFn)
    plt.plot(lNumPoints, np.log(1e-15 + lAcc), label="Gauss-Patterson", marker='s')

    # Gauss-Exotic.
    def pExoGridFn(iDepth):
        exoCT = TasmanianAddons.createExoticQuadratureFromFunction(iDepth + 1, fShift, pWeightFn, iNref, "Exotic Quadrature", iSymmetric)
        return TasmanianSG.makeGlobalGridCustom(iDim, 0, iDepth, 'qptotal', exoCT)
    lNumPoints, lAcc = getData(pExoGridFn, iMaxPoints, fIntegral, pIntegrandFn)
    plt.plot(lNumPoints, np.log(1e-15 + lAcc), label="Gauss-Exotic", marker='*')

    # Show plots
    plt.title(sTitle)
    plt.xlabel("Number of Points")
    plt.ylabel("Log Error")
    plt.legend()
    plt.show()

# Unit Test.
sincFn = lambda x : np.sinc(10 * x / np.pi)
expNegSqrNormFn = lambda xVec : np.exp(-np.linalg.norm(xVec) ** 2)
iMaxPoints = 1000
iDim = 3
fIntegral = 0.32099682841103033 ** iDim
fShift = 1.0
plotAccuracy("Test", iMaxPoints, iDim, fIntegral, expNegSqrNormFn, sincFn, fShift)
