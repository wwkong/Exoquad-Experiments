# Exotic Quadrature Experiments

This repository contains code for replicating a set of quadrature experiments for integrals with signed measures that are bounded below. It contains two folders, and within each folder there is a README file explaining how to run the appropriate code.

A description of each folder is given below.

## BelykinQuad
Contains an implementation of [Beylkin's Gaussian Quadrature method](https://www.sciencedirect.com/science/article/pii/S1063520312001066), which is used as a benchmark. The implementation is built from scratch in Julia.

## ExoQuad
Contains the [TASMANIAN library](https://tasmanian.ornl.gov/) and additional Python code for generating the exotic quadrature summary plots. The results generated in BelykinQuad are copied from the terminal output into the additional Python scripts in this folder.

The C++ implementation of the exotic quadrature functions also resides here (in `Addons/tsgExoticQuadrature.hpp`). 