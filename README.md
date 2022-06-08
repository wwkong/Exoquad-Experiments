Exotic Quadrature Experiments
===============================
A set of quadrature experiments for integrals with signed measures that are bounded below. Written in Julia and Python.  The quadrature method itself can be found in the [TASMANIAN](https://tasmanian.ornl.gov/) library.

Table of Contents
-----------------
* [Requirements](#requirements)
* [Usage](#usage)

Requirements
------------
These experiments require the following to run:
* [CMake](https://cmake.org/) 3.15+
* [Julia](https://julialang.org/) 1.7+
* [Python](https://www.python.org/) 3.8+

Usage
-----
1. Clone this GitHub repository into a source file directory `${SOURCE_DIRECTORY}`, and create a build directory `${BUILD_DIRECTORY}`, e.g., 
```bash
cd ${SOURCE_DIRECTORY}
git clone https://github.com/wwkong/Exoquad-Experiments.git
mkdir ${BUILD_DIRECTORY}
```

2. Create the build files, e.g.,
```bash
cd ${BUILD_DIRECTORY}
cmake -S ${SOURCE_DIRECTORY} -B ${BUILD_DIRECTORY}
make -j
```

3. Generate the data files in the build directory using
```bash
cd ${BUILD_DIRECTORY}
python scripts/GenerateData.py
```
This will also install and compile any necessary Julia packages.

4. Generate the SVG experiment plots using
```bash
cd ${BUILD_DIRECTORY}
python scripts/GeneratePlots_1D.py
python scripts/GeneratePlots_nD1.py
python scripts/GeneratePlots_nD2.py
python scripts/GeneratePlots_nD3.py
```
These plots can be found in `${BUILD_DIRECTORY}/res/`.
