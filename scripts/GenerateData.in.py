import os

# Beylkin experiments (Julia)
os.system("julia --project=\"${CMAKE_CURRENT_BINARY_DIR}/..\" \"${CMAKE_CURRENT_BINARY_DIR}\"/beylkin_antenna.jl")
