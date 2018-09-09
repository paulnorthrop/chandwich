# chandwich 1.0.0.9000

## New features

* `adjust_loglik` has an additional arguments `mle`, `H` and  `V` that allow the user the option to supply the MLE, the Hessian of the independence loglikelihood (H) and the variance of the vector of cluster-specific contributions to the score vector (V), each evaluated at the MLE, rather than estimating these within `adjust_loglik`.

* An anova S3 method for class "chandwich" has been added.  This compares two or more nested models using adjusted likelihood ratio tests of successive pairs of models, using `compare_models()`.

* A confint S3 method for class "chandwich" has been added.  This is based on a fairly trivial call to `conf_inetrvals()`.

* S3 methods `coef` and `vcov` for class "chandwich" have been added.

## Bug fixes and minor improvements

* A bug in `compare_models()` has been fixed.  The bug resulted an error in cases where the argument `larger` corresponded to a simplication of the full model in which element i of the parameter was fixed but some element i+n, for n > 0, was not fixed.

* Estimated unadjusted (VC) and adjusted (adjVC) variance-covariance matrices of the free model parameters are now available in the object returned by `adjust_loglik()`.

* The documentation of the argument `approx` to `compare_models()` has been edited to make it clearer that if `smaller` is not supplied then `approx = FALSE` will be used regardless of any value supplied for `approx` in the call to `compare_models()`.

* If parameter names are supplied to `adjust_loglik()` (via `par_names`) but `fixed_pars` is numeric then the names of the parameters in `fixed_pars` are now also inferred in the case where a larger model is not supplied via `larger`.  This means that the output from `compare_models()` will now use the parameter name, rather than the parameter number.

* If a numeric `fixed_pars` is supplied to `compare_models()` then the names of the parameters in `fixed_pars` are inferred, if they are available in the supplied object `larger`.

* The summary method for class "evpost" is now set up according to Section 8.1 of the R FAQ at (https://cran.r-project.org/doc/FAQ/R-FAQ.html).

* In the Introducing chandwich vignette a typo in the definition of HA has been corrected.  The expression given is for the inverse of HA, not for HA.
