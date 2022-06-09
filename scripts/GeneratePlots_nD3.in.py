import numpy as np
import sys

res_folder = @Exoquad_res_folder@

@Tasmanian_python_import@
import Tasmanian

@Exoquad_src_python_import@
@Exoquad_script_python_import@

# QuadSinc experiments.
import QuadSincExp
for iDim in range(2, 5, 2):
    QuadSincExp.write_plot(20, iDim, 0.15458767050820131, 0.1548241701043440, res_folder)
