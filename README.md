# R Package for Sparse + Low-Rank Inverse Covariance Estimation (SLICE)
This library provides implementation of the Sparse + Low-Rank Inverse Covariance Estimation (SLICE) method for Gaussian graphical model estimation.

## Installation
The package can be installed from this repository with the following code in `R`:

```
# Install devtools if needed, then load
if (!require("devtools", quietly = TRUE))
  install.packages("devtools")
library(devtools)

# Install SLICE
devtools::install_github("mondrus96/slice")
```

## Usage
The main function is `slice` which takes the following form:

```
slice(Sigma, rho, r, Sest = "glasso", tol = 0.001, maxiter = 100)
```

### Arguments

- **Sigma**: A matrix. The input covariance matrix.
- **rho**: A numeric. Regularization parameter for sparse estimator.
- **r**: An integer. Rank for latent component.
- **Sest**: A character string. Type of sparse estimator to use, default = `"glasso"` [(Friedman et al, 2008)](https://academic.oup.com/biostatistics/article/9/3/432/224260). Other choices include `"gscad"` [(Fan et al., 2009)](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-3/issue-2/Network-exploration-via-the-adaptive-LASSO-andSCAD-penalties/10.1214/08-AOAS215.full), `"clime"` [(Cai et al., 2011)](https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.tm10155), and `"huge_glasso"` [(Zhao et al., 2012)](https://jmlr.org/beta/papers/v13/zhao12a.html).
- **tol**: A numeric. Tolerance for algorithm, default = 1e-3.
- **maxiter**: An integer. Maximum number of iterations, default = 100.

### Returns

An S3 class `slice` object with:

- **S**: A matrix corresponding to the estimated sparse component.
- **L**: A matrix corresponding to the estimated latent component.
- **rho**: A numeric of the regularization parameter used for the sparse component.
- **r**: An integer of the rank used for the latent component.
- **misc**: Contains additional outputs related to the convergence of the algorithm.

## Examples
We provide an example of `slice` a given covariance matrix, `Sigma`, sparsity tuning parameter `rho`, and rank of latent matrix `r`. Using the default `glasso` method:

```
out <- slice(Sigma, rho, r) # Run SLICE
```

Alternate penalties of the sparse component can be specified. For example, for SCAD:

```
out <- slice(Sigma, rho, r, Sest = "gscad") # Run SLICE with SCAD
```