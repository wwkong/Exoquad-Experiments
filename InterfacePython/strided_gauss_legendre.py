import numpy as np
import math

import Tasmanian
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pandas as pd

def getGaussLegendreSubset(dim, depth, max_level, stride):
# Create a grid containing a subset of the Gauss-Legendre rules.
    grid = Tasmanian.TasmanianSparseGrid()
    numLevels = 0
    numNodes = np.array([], dtype=np.int64)
    precision = np.array([], dtype=np.int64)
    weights = []
    nodes = []
    for level in range(0, max_level, stride):
        numLevels += 1
        grid.makeGlobalGrid(1, 0, level, "level", "gauss-legendre")
        numNodes = np.append(numNodes, grid.getNumPoints())
        precision = np.append(precision, 2 * grid.getNumPoints() - 1)
        weights.append(grid.getQuadratureWeights())
        nodes.append(grid.getPoints().flatten())
    # Create the subset CustomTabulated instance and the derived global grid
    subCT = Tasmanian.makeCustomTabulatedFromData(numLevels, numNodes, precision, nodes, weights, "Gauss-Legendre Subset")
    grid.makeGlobalGridCustom(dim, 0, depth, "qptotal", subCT)
    return(grid)

def testSubsetAccuracy(dim, depth, max_level, stride):
# Test the accuracy of a subset grid for the function:
#   (x1,...,xd) -> exp(-|x|^2) * cos(x1) * ... * cos(xd).
    grid = getGaussLegendreSubset(dim, depth, max_level, stride)
    weights = grid.getQuadratureWeights()
    nodes = grid.getPoints()
    approx_integral = 0.0
    fn_vals = np.exp(-np.linalg.norm(nodes - 0.2, axis=1) ** 2)
    approx_integral = np.dot(fn_vals, weights)
    exact_integral = pow(1.464414614451514, dim)
    log10_err = math.log10(1e-15 + abs(approx_integral - exact_integral))
    return(dim, depth, stride, nodes.size, log10_err)

def plot_errors(dim, strides, max_obs=3e5, log10_err_tol=-12):
    # Plots the log10 error across different depths and strides, for a fixed dim
    data_array = np.empty((0, 5))
    for stride in strides:
        log10_err = float("inf")
        obs = 0
        depth = 0
        while (log10_err > log10_err_tol and obs < max_obs):
            print([dim, depth, stride])
            row = np.asarray([testSubsetAccuracy(dim, depth, depth * stride + 1, stride)])
            data_array = np.append(data_array, row, axis=0)
            depth += stride
            obs = row[0, 3]
            log10_err = row[0, 4]
    df = pd.DataFrame(dict(dim=data_array[:,0], depth=data_array[:,1], stride=data_array[:,2].astype(int),
                           numPoints=data_array[:,3], log10_err=data_array[:,4]))
    groups = df.groupby('stride')
    plt.rc('font', size=12)
    fig, ax = plt.subplots()
    for name, group in groups:
        ax.plot(group.numPoints, group.log10_err, label=name, linestyle='solid', linewidth=0.8, marker='.')
    plt.legend(title='Stride', loc='upper right')
    plt.xlabel('Number of Points')
    plt.ylabel(r'$\log($1E-15 $+ |{\rm Error}|)$')
    plt.title("Log error of the Gauss-Legendre rule for \n $\int_{[-1,1]^" + str(dim) +
              "} \exp(-\|x - (0.2, \ldots, 0.2)\|^2) dx$")
    plt.savefig("gl_dim" + str(dim) + ".svg")
    # plt.show()

# Plot a specific instance
for dim in range(1, 6+1):
    if dim <= 4:
        strides = range(1, 8+1)
    else:
        strides = range(1, 6+1)
    plot_errors(dim, strides)
