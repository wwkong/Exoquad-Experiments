# One-dimensional experiments for the exotic quadrature paper.
import ExoquadUtils
import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

# f(x) = exp(-x^2), rho(x) = sin(a * pi * |x|)
def write_plot(a, fBaseIntegral, dname):
    pIntegrandFn = lambda x : np.exp(-x ** 2)
    pWeightFn = lambda x : np.sin(a * np.pi * np.abs(x))
    sTitle = r'Relative Errors for Frequency $a=' + str(int(a)) +  r'$'
    sFname = dname + '/abs_sin_a' + str(a) + '.svg'
    print('Creating ' + sFname + '...')
    ExoquadUtils.plotAccuracy(sTitle, 25, 1e-14, 1, fBaseIntegral, pIntegrandFn, pWeightFn, 1.0, sFname=sFname, iLocalPolynomial=True)
