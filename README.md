
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Travis-CI Build Status](https://travis-ci.org/paulnorthrop/chandwich.svg?branch=master)](https://travis-ci.org/paulnorthrop/chandwich) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/paulnorthrop/chandwich?branch=master&svg=true)](https://ci.appveyor.com/project/paulnorthrop/chandwich) [![Coverage Status](https://codecov.io/github/paulnorthrop/chandwich/coverage.svg?branch=master)](https://codecov.io/github/paulnorthrop/chandwich?branch=master)

chandwich: Chandler-Bate Sandwich Loglikelihood Adjustment
----------------------------------------------------------

### What does chandwich do?

The `chandwich` package performs adjustments of an independence loglikelihood using a robust sandwich estimator of the parameter covariance matrix, based on the methodology in [Chandler and Bate (2007)](http://dx.doi.org/10.1093/biomet/asm015). This can be used for cluster correlated data when interest lies in the parameters of the marginal distributions.

### A simple example

The main function in the threshr package is `adjust_loglik`. It finds the maximum likelihood estimate (MLE) of model parameters based on an independence loglikelihood in which cluster dependence in the data is ignored. The independence loglikelihood is adjusted in a way that ensures that the Hessian of the adjusted loglikelihood conicides with a robust sandwich estimate of the parameter covariance at the MLE. Three adjustments are available: one in which the independence loglikelihood itself is scaled (vertical scaling) and two others where the scaling is in the parameter vector (horizontal scaling).

``` r
binom_loglik <- function(prob, data) {
  if (prob < 0 || prob > 1) {
    return(-Inf)
  }
  return(dbinom(data[, "y"], data[, "n"], prob, log = TRUE))
}
rat_res <- adjust_loglik(loglik = binom_loglik, data = rats)
```

### Installation

To install this development version from Github use:

``` r
install.packages("devtools")
devtools::install_github("paulnorthrop/chandwich")
```

### Vignette

See `vignette("threshr-vignette", package = "threshr")` for an overview of the package.
