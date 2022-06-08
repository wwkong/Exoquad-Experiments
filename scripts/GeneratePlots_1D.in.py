import sys
import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt

res_folder = @Exoquad_res_folder@

@Tasmanian_python_import@
import Tasmanian

@Exoquad_src_python_import@
@Exoquad_script_python_import@
import ExoquadUtils

# Intro plot
iMaxPoints = 30
fMinErr = 1e-13
iFontSize = 15
fName = res_folder + "/sinc_plot.svg"
print('Creating ' + fName + "...")
plt.rc('axes', titlesize=iFontSize)
plt.rc('axes', labelsize=iFontSize)
plt.rc('xtick', labelsize=iFontSize-2)
plt.rc('ytick', labelsize=iFontSize-2)
plt.rc('figure', titlesize=iFontSize+2)
plt.figure()
for (a, m) in [(1,'o'), (5,'v'), (10,'v'), (25,'D')]:
    fIntegral = 2 * sp.sici(a)[0] / a
    pIntegrandFn = lambda x : np.sinc(a * x / np.pi)
    fErrLower = 1e-20
    pGridFn = lambda iDepth : Tasmanian.makeGlobalGrid(1, 0, iDepth, 'qptotal', 'gauss-legendre')
    lNumPoints, lErr, _ = ExoquadUtils.getData(pGridFn, iMaxPoints, fMinErr, fIntegral, pIntegrandFn)
    plt.plot(lNumPoints, fErrLower + lErr, label=r'$a=' + str(a) + r'$', marker=m)
plt.title("Relative Errors for Different Frequencies")
plt.xlabel("Number of Points")
plt.ylabel("Relative Error")
plt.yscale("log")
plt.grid(linestyle="--", linewidth=0.25)
plt.legend()
plt.savefig(fName)

# 1D Antenna experiments (a=10,20).
import AntennaExp
antenna_a10_fname = res_folder + "/beylkin_antenna_data_a10.txt"
with open(antenna_a10_fname, "r") as f:
    antenna_a10_data = np.array([[float(num) for num in line.split('\t')] for line in f])
AntennaExp.write_plot(10, antenna_a10_data, res_folder)
antenna_a20_fname = res_folder + "/beylkin_antenna_data_a20.txt"
with open(antenna_a20_fname, "r") as f:
    antenna_a20_data = np.array([[float(num) for num in line.split('\t')] for line in f])
AntennaExp.write_plot(20, antenna_a20_data, res_folder)

# 1D AbsSin experiments (a=2,5).
import AbsSinExp
AbsSinExp.write_plot(2, 0.22869582997457813, res_folder)
AbsSinExp.write_plot(5, 0.17482639269260326, res_folder)

