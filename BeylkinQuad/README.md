# Beylkin's Gaussian Quadrature
Implementation of methods by Beylkin et al. on Generalized Gaussian Quadrature for Bandlimited Exponentials

To run the experiments, use the following command in this folder:
```console
julia --project=. test/runtests.jl
```
This will print a matrix to the terminal. The first column of the matrix is the number of points used in the quadrature, while the second column is the relative error given by

```python
relative_error = abs(approx_integral - true_integral) / abs(true_integral)
```