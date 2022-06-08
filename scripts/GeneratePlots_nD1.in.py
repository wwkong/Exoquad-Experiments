import numpy as np
import sys

res_folder = @Exoquad_res_folder@

@Tasmanian_python_import@
import Tasmanian

@Exoquad_src_python_import@
@Exoquad_script_python_import@

# SqrSinc experiments.
import SqrSincExp
for iDim in range(2, 5, 2):
    SqrSincExp.write_plot(10, iDim, 0.32099682841103033, res_folder)
