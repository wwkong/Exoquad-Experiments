# One-dimensional experiments for the exotic quadrature paper.
import ExoquadUtils
import Tasmanian
import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

# f(x) = exp(-x^2), rho(x) = sinc(a * pi)
def write_plot(a, fBaseIntegral, dname):
    sFname = dname + '/sinc_a' + str(a) + '.svg'
    print('Creating ' + sFname + '...', end='', flush=True)
    pIntegrandFn = lambda x : np.exp(-x ** 2)
    pWeightFn = lambda x : np.sinc((a * x) / np.pi)
    sTitle = r'Relative Errors for Frequency $a=' + str(int(a)) +  r'$'
    ExoquadUtils.plotAccuracy(sTitle, 30, 1E-12, 1, fBaseIntegral, pIntegrandFn, pWeightFn, 1.0, sFname=sFname, iNref=1000)
    print("DONE")
