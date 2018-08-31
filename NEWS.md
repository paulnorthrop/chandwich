# chandwich 1.0.0.9000

## New features

* `adjust_loglik` has an additional arguments `mle`, `H` and  `V` that allow the user the option to supply the MLE, the Hessian of the independence loglikelihood (H) and the variance of the vector of cluster-specific contributions to the score vector (V), each evaluated at the MLE, rather than estimating these within `adjust_loglik`.

## Bug fixes and minor improvements

The documentation of the argument `approx` to `compare_models()` has been edited to make it clearer that if `smaller` is not supplied then `approx = FALSE` will be used regardless of any value supplied for `approx` in the call to `compare_models()`.

In the Introducing chandwich vignette a typo in the definition of HA has been corrected.  The expression given is for the inverse of HA, not for HA.
