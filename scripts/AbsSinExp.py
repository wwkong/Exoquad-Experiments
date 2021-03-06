# One-dimensional experiments for the exotic quadrature paper.
import ExoquadUtils
import Tasmanian
import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

# f(x) = exp(-x^2), rho(x) = sin(a * pi * |x|)
def write_plot(a, fBaseIntegral, dname):
    sFname = dname + '/abs_sin_a' + str(a) + '.svg'
    print('Creating ' + sFname + '...', end='', flush=True)
    pIntegrandFn = lambda x : np.exp(-x ** 2)
    pWeightFn = lambda x : np.sin(a * np.pi * np.abs(x))
    localPolyGrid = ExoquadUtils.getRefinedSurrogate(pWeightFn, 7, 'fds', 1E-7)
    sTitle = r'Relative Errors for Frequency $a=' + str(int(a)) +  r'$'
    ExoquadUtils.plotAccuracy(sTitle, 60, 1E-7, 1, fBaseIntegral, pIntegrandFn, pWeightFn, 1.0, sFname=sFname,
                              iLocalPolynomial=True, oRefGrid=localPolyGrid)
    print("DONE")
