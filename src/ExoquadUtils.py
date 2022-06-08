# One-dimensional experiments for the exotic quadrature paper.
import Tasmanian
import TasmanianAddons
import numpy as np
import matplotlib.pyplot as plt
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

def plotAccuracy(sTitle, iMaxPoints, fMinErr, iDim, fIntegral, pIntegrandFn, pWeightFnReal, fShift, pWeightFnIm=[],
                 lGBNumPoints=np.array([]), lGBRelativeErr=np.array([]), sFname="plot.svg", iSymmetric=False, iNref=500,
                 iLocalPolynomial=False):
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
    fig = plt.figure()
    if pWeightFnIm:
        pCombinedFn = lambda x : pIntegrandFn(x) * np.prod(pWeightFnReal(x) +  pWeightFnIm(x) * 1j)
    else:
        pCombinedFn = lambda xVec : pIntegrandFn(xVec) * np.prod(pWeightFnReal(xVec))
    fErrLower = 1e-20

    # Wrappers.
    def pGridWrapper(iDepth, sType):
        return Tasmanian.makeGlobalGrid(iDim, 0, iDepth, 'qptotal', sType)

    # Gauss-Legendre (all).
    pGridFn = lambda iDepth : pGridWrapper(iDepth, 'gauss-legendre')
    lNumPoints, lErr, _ = getData(pGridFn, iMaxPoints, fMinErr, fIntegral, pCombinedFn)
    plt.plot(lNumPoints, fErrLower + lErr, label="GL", marker='o')

    # Gauss-Legendre (odd).
    pGridFn = lambda iDepth : pGridWrapper(iDepth, 'gauss-legendre-odd')
    lNumPoints, lErr, _ = getData(pGridFn, iMaxPoints, fMinErr, fIntegral, pCombinedFn)
    plt.plot(lNumPoints, fErrLower + lErr, label="GL-odd", marker='v')

    # Gauss-Patterson.
    pGridFn = lambda iDepth : pGridWrapper(iDepth, 'gauss-patterson')
    lNumPoints, lErr, _ = getData(pGridFn, iMaxPoints, fMinErr, fIntegral, pCombinedFn)
    plt.plot(lNumPoints, fErrLower + lErr, label="GP", marker='s')

    # Gauss-Exotic.
    def pExoGridFn(iDepth, pWeightFn):
        exoCT = TasmanianAddons.createExoticQuadratureFromFunction(iDepth + 1, fShift, pWeightFn, iNref, "Exotic Quadrature", iSymmetric)
        return Tasmanian.makeGlobalGridCustom(iDim, 0, iDepth, 'qptotal', exoCT)
    if pWeightFnIm:
        pExoGridFnReal = lambda iDepth : pExoGridFn(iDepth, pWeightFnReal)
        pExoGridFnIm = lambda iDepth : pExoGridFn(iDepth, pWeightFnIm)
        dummy = 1.0
        lNumPoints, _, lEstReal = getData(pExoGridFnReal, iMaxPoints, fErrLower, dummy, pIntegrandFn)
        _, _, lEstIm = getData(pExoGridFnIm, iMaxPoints, fErrLower, dummy, pIntegrandFn)
        lErr = abs(lEstReal + lEstIm * 1j - fIntegral) / abs(fIntegral)
        # Filter
        lNumPoints = lNumPoints[lErr >= fMinErr]
        lErr = lErr[lErr >= fMinErr]
    else:
        pExoGridFnReal = lambda iDepth : pExoGridFn(iDepth, pWeightFnReal)
        lNumPoints, lErr, _ = getData(pExoGridFnReal, iMaxPoints, fMinErr, fIntegral, pIntegrandFn)
    plt.plot(lNumPoints, fErrLower + lErr, label="Exotic", marker='*')

    # Local Polynomial
    if iLocalPolynomial:
        pGridFn = lambda iDepth : Tasmanian.makeLocalPolynomialGrid(iDim, 0, iDepth, iOrder=2, sRule="localp")
        lNumPoints, lErr, _ = getData(pGridFn, iMaxPoints, fMinErr, fIntegral, pCombinedFn)
        plt.plot(lNumPoints, fErrLower + lErr, label="LP", marker='+')

    # Gauss-Bandlimited
    if lGBNumPoints.any():
        plt.plot(lGBNumPoints, fErrLower + lGBRelativeErr, label="GB", marker='D')

    # Show plots.
    plt.title(sTitle)
    plt.xlabel("Number of Points")
    plt.ylabel("Relative Error")
    plt.yscale("log")
    plt.grid(linestyle="--", linewidth=0.25)
    if lGBNumPoints.any() or iLocalPolynomial:
        # Re-arrange legend.
        handles, labels = plt.gca().get_legend_handles_labels()
        order = [4,0,1,2,3]
        plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc="upper right")
    else:
        plt.legend(loc="upper right")
    plt.savefig(sFname)
