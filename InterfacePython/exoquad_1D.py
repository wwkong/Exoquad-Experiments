# One-dimensional experiments for the exotic quadrature paper.
import ExoquadUtils
import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

# ================================================================================================================================
# Experiment #1, f(x) = exp(im * a * x), rho(x) = I0(pi * sqrt(1 - x^2)) / 2
# ================================================================================================================================
# Raw arrays from Julia (a=10).
lGBNumPoints = [8,9,10,12,13,14,15,16,17,18,19,20,21,22]
lGBRelativeErrs = [41.0496,56.369,57.0769,3.87772,0.587056,0.0466392,0.023792,0.0265805,0.02104,0.0176328,0.0130911,0.00341601,
                    0.00912594,0.0602533]
# Python arrays (a=10).
a = 10
pIntegrandFn = lambda x : sp.iv(0, np.pi * np.sqrt(1 - x ** 2)) / 2
pWeightFnReal = lambda x : np.cos(a * np.pi * x)
pWeightFnIm = lambda x : np.sin(a * np.pi * x)
fBaseIntegral = np.sinc(np.sqrt((a * np.pi) ** 2 - np.pi ** 2) / np.pi)
sTitle = r'Relative Errors for Frequency $a=' + str(int(a)) +  r'$'
sFname = 'antenna_a10.svg'
print('Creating ' + sFname + '...')
ExoquadUtils.plotAccuracy(sTitle, 25, 1e-14, 1, fBaseIntegral, pIntegrandFn, pWeightFnReal, 1.0,
                           pWeightFnIm=pWeightFnIm, lGBNumPoints=lGBNumPoints, lGBRelativeErr=lGBRelativeErrs, sFname=sFname)

# Raw arrays from Julia (a=20).
lGBPoints = [17,18,19,20,22,23,24,25,27,27,28,29,30,31,32,33,34]
lGBErrs = [85.7029,75.4862,1.32635,137.758,16.786,3.58022,0.441072,0.148161,0.181965,0.181965,0.172947,0.129785,0.111163,
           0.105705,0.0516552,0.0990842,0.796096]
# Python arrays (a=20).
a = 20
pIntegrandFn = lambda x : sp.iv(0, np.pi * np.sqrt(1 - x ** 2)) / 2
pWeightFnReal = lambda x : np.cos(a * np.pi * x)
pWeightFnIm = lambda x : np.sin(a * np.pi * x)
fBaseIntegral = np.sinc(np.sqrt((a * np.pi) ** 2 - np.pi ** 2) / np.pi)
sTitle = r'Relative Errors for Frequency $a=' + str(int(a)) +  r'$'
sFname = 'antenna_a20.svg'
print('Creating ' + sFname + '...')
ExoquadUtils.plotAccuracy(sTitle, 25, 1e-14, 1, fBaseIntegral, pIntegrandFn, pWeightFnReal, 1.0,
                          pWeightFnIm=pWeightFnIm, lGBNumPoints=lGBNumPoints, lGBRelativeErr=lGBRelativeErrs, sFname=sFname)

# ================================================================================================================================
# Experiment #2, f(x) = exp(-x^2), rho(x) = sin(a * pi * |x|)
# ================================================================================================================================

a = 2
pIntegrandFn = lambda x : np.exp(-x ** 2)
pWeightFn = lambda x : np.sin(a * np.pi * np.abs(x))
fBaseIntegral = 0.22869582997457813
sTitle = r'Relative Errors for Frequency $a=' + str(int(a)) +  r'$'
sFname = 'abs_sin_a' + str(a) + '.svg'
print('Creating ' + sFname + '...')
ExoquadUtils.plotAccuracy(sTitle, 25, 1e-14, 1, fBaseIntegral, pIntegrandFn, pWeightFn, 1.0, sFname=sFname, iLocalPolynomial=True)

a = 5
pIntegrandFn = lambda x : np.exp(-x ** 2)
pWeightFn = lambda x : np.sin(a * np.pi * np.abs(x))
fBaseIntegral = 0.17482639269260326
sTitle = r'Relative Errors for Frequency $a=' + str(int(a)) +  r'$'
sFname = 'abs_sin_a' + str(a) + '.svg'
print('Creating ' + sFname + '...')
ExoquadUtils.plotAccuracy(sTitle, 25, 1e-14, 1, fBaseIntegral, pIntegrandFn, pWeightFn, 1.0, sFname=sFname, iLocalPolynomial=True)

