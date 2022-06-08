import numpy as np
import sys

res_folder = @Exoquad_res_folder@

@Tasmanian_python_import@
import Tasmanian

@Exoquad_src_python_import@
@Exoquad_script_python_import@

# SumSinc experiments.
import SumSincExp
for iDim in range(2, 5, 2):
    SumSincExp.write_plot(10, iDim, 0.34005259521041936, res_folder)
