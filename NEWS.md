# chandwich 1.1.4.9000

## Bug fixes and minor improvements

* In the vignette references section DOIs are rendered with URLs hyperlinked.

# chandwich 1.1.4

## Bug fixes and minor improvements

* The function `conf_intervals()` now stops evaluating the profile loglikelihood when the required confidence limits have been estimated.

* In `conf_interals()` the initial estimates used in the profile loglikelihood were not reset in between the searches for the lower confidence limits and the upper limits.  In some cases this may cause problems for the optimisations.  This has been corrected: initial estimates are now reset. 

* In `confint.chandwich()` the argument `profile` can be used to choose whether to return confidence intervals based on an (adjusted) profile loglikelihood or based on approximate large sample normal theory.

* The `anova.chandwich()` function has been modified to enable the comparison of nested (adjusted) model objects, where the nesting has not been explicitly created using the argument `fixed_pars` to `adjust_loglikelihood()`.

* The documentation of `conf_intervals()` has been updated: the general nature of the argument `num` is explained and in **Details** advice is given on what to do if one or more of the confidence limits are not found using the default arguments.

* In the description of the argument `parm` in the documentation of `confint.chandwich()`, `which_pars` has been corrected to `parm` twice.

# chandwich 1.1.3

## Bug fixes and minor improvements

* The object returned from `adjust_loglik()` has an extra attribute `loglikVecMLE`, which is a vector containing the contributions of individual observations to the independence log-likelihood evaluated at the MLE.

* `isTRUE()` is used in the GEV example to avoid potential problems owing to the possibility of NAs.

# chandwich 1.1.2

## Bug fixes and minor improvements

* A typo meant that the text in the Value section of the description of `confint.chandwich()` was cut short.  This has been corrected.

* In `adjust_loglik()` if `p = 1` and the user supplies a scalar `H` instead of a matrix then an error is thrown by `dimnames()` just before the results are returned.  The calculations are correct, so the code has been modified trivially to avoid the error.

* In the function returned from `adjust_loglik` the code to perform the horizontal adjustment of the loglikelihood has been changed to `x_star <- mle + as.vector(C %*% (x - mle))`: using `as.vector()` avoids a potential warning by ensuring vector + vector, not vector + matrix.

* There were bugs in `plot.chandwich()`, `plot.confreg()` and `plot.confint()` that produced an error if the user tried to set their own legend using the `legend` argument.  These bugs have been corrected.

# chandwich 1.1.1

## Bug fixes and minor improvements

* The attribute `nobs` has been added to the object returned from `adjust_loglik()` and as an attribute to the object returned from `logLik.chandwich()`. 

* In `compare_models()` the parameter names (if any) are passed to the (adjusted) loglikelihood function, in case they are required inside the loglikelihood function.

* There was a bug in the plot method for objects of class "confreg" returned from conf_region(): if the parameters had not been named by the user then ? appeared twice in the console, requiring the user to press return twice before the plot as produced.  This has been corrected. 

# chandwich 1.1.0

## New features

* `adjust_loglik` has an additional arguments `mle`, `H` and  `V` that allow the user the option to supply the MLE, the Hessian of the independence loglikelihood (H) and the variance of the vector of cluster-specific contributions to the score vector (V), each evaluated at the MLE, rather than estimating these within `adjust_loglik`.

* An anova S3 method for class "chandwich" has been added.  This compares two or more nested models using adjusted likelihood ratio tests of successive pairs of models, using `compare_models()`.

* A confint S3 method for class "chandwich" has been added.  This is based on a fairly trivial call to `conf_intervals()`.

* S3 methods `coef`, `vcov` and `logLik` for class "chandwich" have been added.

## Bug fixes and minor improvements

* A bug in `compare_models()` has been fixed.  The bug resulted an error in cases where the argument `larger` corresponded to a simplification of the full model in which element i of the parameter was fixed but some element i+n, for n > 0, was not fixed.

* Estimated unadjusted (VC) and adjusted (adjVC) variance-covariance matrices of the free model parameters are now available in the object returned by `adjust_loglik()`.

* The documentation of the argument `approx` to `compare_models()` has been edited to make it clearer that if `smaller` is not supplied then `approx = FALSE` will be used regardless of any value supplied for `approx` in the call to `compare_models()`.

* If parameter names are supplied to `adjust_loglik()` (via `par_names`) but `fixed_pars` is numeric then the names of the parameters in `fixed_pars` are now also inferred in the case where a larger model is not supplied via `larger`.  This means that the output from `compare_models()` will now use the parameter name, rather than the parameter number.

* If a numeric `fixed_pars` is supplied to `compare_models()` then the names of the parameters in `fixed_pars` are inferred, if they are available in the supplied object `larger`.

* The summary method for class "evpost" is now set up according to Section 8.1 of the R FAQ at (https://cran.r-project.org/doc/FAQ/R-FAQ.html).

* In the Introducing chandwich vignette a typo in the definition of HA has been corrected.  The expression given is for the inverse of HA, not for HA.
