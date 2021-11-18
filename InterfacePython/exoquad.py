# Experiments for the exotic quadrature paper.
import TasmanianSG
import TasmanianAddons
import numpy as np
import math

def createGrids(iDim, iDepth, oWeightFn, fShift, iSymmetric = False, iNref = 100):
# Create Gauss-Legendre, Gauss-Legendre (odd), Gauss-Patterson, and Exotic Quadrature grids, and output them as a list of
# TasmanianSparseGrid objects.
    glGrid = TasmanianSG.makeGlobalGrid(iDim, 0, iDepth, 'qptotal', 'gauss-legendre')
    glOddGrid = TasmanianSG.makeGlobalGrid(iDim, 0, iDepth, 'qptotal', 'gauss-legendre-odd')
    gpGrid = TasmanianSG.makeGlobalGrid(iDim, 0, iDepth, 'qptotal', 'gauss-patterson')
    exoCT = TasmanianAddons.createExoticQuadratureFromFunction(iDepth, fShift, oWeightFn, iNref, "Exotic Quadrature", iSymmetric)
    exoGrid = TasmanianSG.makeGlobalGridCusto(iDim, 0, iDepth, 'qptotal', exoCT)
    return [glGrid, glOddGrid, gpGrid, exoGrid]

def plotGrids(lGrids, lNames, fIntegral, oIntegrandFn, oWeightFn):
# Plots the accuracy of a list of grids against a reference integral and a given a 1D weight function, nD integrand function.
    return

