# One-dimensional experiments for the exotic quadrature paper.
import ExoquadUtils
import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

# f(x) = exp(im * a * x), rho(x) = I0(pi * sqrt(1 - x^2)) / 2
def write_plot(a, beylkin_data, dname):
    pIntegrandFn = lambda x : sp.iv(0, np.pi * np.sqrt(1 - x ** 2)) / 2
    pWeightFnReal = lambda x : np.cos(a * np.pi * x)
    pWeightFnIm = lambda x : np.sin(a * np.pi * x)
    fBaseIntegral = np.sinc(np.sqrt((a * np.pi) ** 2 - np.pi ** 2) / np.pi)
    sTitle = r'Relative Errors for Frequency $a=' + str(int(a)) +  r'$'
    sFname = dname + "/antenna_a" + str(int(a)) + ".svg"
    print('Creating ' + sFname + '...')
    ExoquadUtils.plotAccuracy(sTitle, 25, 1e-8, 1, fBaseIntegral, pIntegrandFn, pWeightFnReal, 1.0,
                              pWeightFnIm=pWeightFnIm, lGBNumPoints=beylkin_data[:,0], lGBRelativeErr=beylkin_data[:,1],
                              sFname=sFname, iNref=2000)
