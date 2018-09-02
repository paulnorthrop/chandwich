# chandwich 1.0.0.9000

## New features

* `adjust_loglik` has an additional arguments `mle`, `H` and  `V` that allow the user the option to supply the MLE, the Hessian of the independence loglikelihood (H) and the variance of the vector of cluster-specific contributions to the score vector (V), each evaluated at the MLE, rather than estimating these within `adjust_loglik`.

## Bug fixes and minor improvements

* The documentation of the argument `approx` to `compare_models()` has been edited to make it clearer that if `smaller` is not supplied then `approx = FALSE` will be used regardless of any value supplied for `approx` in the call to `compare_models()`.

* If parameter names are supplied to `adjust_loglik()` (via `par_names`) but `fixed_pars` is numeric then the names of the parameters in `fixed_pars` are now also inferred in the case where a larger model is not supplied via `larger`.  This means that the output from `compare_models()` will now use the parameter name, rather than the parameter number.

* If a numeric `fixed_pars` is supplied to `compare_models()` then the names of the parameters in `fixed_pars` are inferred, if they are available in the supplied object `larger`.

* The summary method for class "evpost" is now set up according to Section 8.1 of the R FAQ at (https://cran.r-project.org/doc/FAQ/R-FAQ.html).

* In the Introducing chandwich vignette a typo in the definition of HA has been corrected.  The expression given is for the inverse of HA, not for HA.
