# Exotic Quadrature
Additional Python code for running exotic quadrature experiments in the TASMANIAN library

To run the experiments:

1. Build the library using CMake with Python support enabled, e.g.,
```bash
#!/bin/sh

BUILD_DIR=/tmp/ramdisk/TASMANIAN_Build
INSTALL_DIR=/tmp/ramdisk/TASMANIAN_Install
SOURCE_DIR=~/Projects/TASMANIAN

if [ -d "$BUILD_DIR" ]; then rm -rf $BUILD_DIR; fi
mkdir $BUILD_DIR

if [ -d "$INSTALL_DIR" ]; then rm -Rf $INSTALL_DIR; fi
mkdir $INSTALL_DIR

cd $BUILD_DIR
cmake -D CMAKE_INSTALL_PREFIX:PATH="$INSTALL_DIR" \
      -D CMAKE_EXPORT_COMPILE_COMMANDS=1 \
      -D CMAKE_CXX_COMPILER="/bin/clang++" \
      -D CMAKE_CXX_FLAGS="-O3 -Wall -Wextra -Wshadow" \
      -D Tasmanian_ENABLE_RECOMMENDED=ON \
      $SOURCE_DIR
make -j
```

2. In the build directory, generate the Python plots with the following commands:
```console
python $BUILD_DIR/Python/exoquad_1D.py
python $BUILD_DIR/Python/exoquad_nD.py
```
Specifically, this will generate a set of .svg files containing the plots.