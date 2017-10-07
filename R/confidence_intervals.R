# ============================== profile_loglik ===============================

#' Profile loglikelihood
#'
#' Calculates the profile loglikelihood for a subset of the model parameters.
#'
#' @param object An object of class \code{"chandwich"} returned by
#'   \code{adjust_loglik}.
#' @param prof_pars A vector specifying the subset of the (unfixed) parameters
#'   over which to profile. Can be either a numeric vector, specifying indices
#'   of the components of the \strong{full} parameter vector, or a character
#'   vector of parameter names, which must be a subset of those supplied in
#'   \code{par_names} in the call to \code{\link{adjust_loglik}} that produced
#'   \code{object}.
#'
#'   \code{prof_pars} must not have any parameters in common with
#'   \code{attr(object, "fixed_pars")}.  \code{prof_pars} must not contain
#'   all of the unfixed parameters, i.e. there is no point in profiling over
#'   all the unfixed parameters.
#' @param prof_vals A numeric vector.  Values of the parameters in
#'   \code{prof_pars}.  If \code{prof_vals = NULL} then the MLEs stored
#'   in \code{object} of the parameters \code{prof_pars} are used.
#' @param init A numeric vector of initial estimates of the values of the
#'   parameters that are not fixed and are not in \code{prof_pars}.
#'   Should have length \code{attr(object, "p_current") - length(prof_pars)}.
#'   If \code{init} is \code{NULL} or is of the wrong length then the
#'   relevant components from the MLE stored in \code{object} are used.
#' @param type A character scalar.  The argument \code{type} to the function
#'   returned by \code{\link{adjust_loglik}}, that is, the type of adjustment
#'   made to the independence loglikelihood function.
#' @param ... Further arguments to be passed to \code{\link[stats]{optim}}.
#'   These may include \code{gr}, \code{method}, \code{lower}, \code{upper}
#'   or \code{control}.
#' @return A numeric vector of length 1.  The value of the profile
#'   loglikelihood.
#' @examples
#' profile_loglik(large, prof_pars = "xi1", prof_vals = -0.1)
profile_loglik <- function(object, prof_pars = NULL, prof_vals = NULL,
                           init = NULL, type = c("vertical", "cholesky",
                                                 "dilation", "none"), ...) {
  type <- match.arg(type)
  if (is.null(prof_pars)) {
    stop("prof_pars must be supplied")
  }
  fixed_pars <- attr(object, "fixed_pars")
  fixed_at <- attr(object, "fixed_at")
  full_par_names <- attr(object, "full_par_names")
  if (!all(prof_pars %in% full_par_names)) {
    stop("prof_pars is not a subset of ", deparse(full_par_names))
  }
  if (is.character(prof_pars)) {
    if (is.null(full_par_names)) {
      stop("prof_pars can be character only if par_names is supplied")
    }
    temp <- prof_pars
    prof_pars <- which(full_par_names %in% prof_pars)
    names(prof_pars) <- temp
  } else {
    if (!is.null(full_par_names)) {
      names(prof_pars) <- full_par_names[prof_pars]
    }
  }
  rest_pars <- (1:attr(object, "p_full"))[-c(fixed_pars, prof_pars)]
  if (!is.null(full_par_names)) {
    names(rest_pars) <- full_par_names[rest_pars]
  }
  if (any(prof_pars %in% fixed_pars)) {
    stop("prof_pars and attr(object,''fixed_pars'') have parameters in common")
  }
  if (length(rest_pars) == 0) {
    stop("There is no point in profiling over all the unfixed parameters.")
  }
  # The MLE (including any fixed parameters) of the full parameter vector
  full_mle <- attr(object, "res_MLE")
  if (is.null(prof_vals)) {
    prof_vals <- full_mle[prof_pars]
  }
  p_r <- length(rest_pars)
  # Initial estimate: the MLE with fixed_pars set at the values in fixed_at
  if (is.null(init) || length(init) != p_r) {
    init <- full_mle[rest_pars]
  }
  # Extract arguments to be passed to optim()
  optim_args <- list(...)
  # Use "BFGS", unless the user has chosen the method or if p_r = 1 and they
  # have (inappropriately) chosen "Nelder-Mead" when p_r = 1
  if (is.null(optim_args$method)) {
    optim_args$method <- "BFGS"
  } else if (p_r == 1 & optim_args$method == "Nelder-Mead") {
    optim_args$method <- "BFGS"
  }
  # Set up a function to perform the profiling
  # Fixed_pars are dealt with inside the function returned by adjust_loglik()
  p <- attr(object, "p_current")
  neg_prof_loglik <- function(x) {
    pars <- numeric(p)
    pars[prof_pars] <- prof_vals
    pars[rest_pars] <- x
    return(-do.call(object, list(pars, type = type)))
  }
  # L-BFGS-B and Brent don't like Inf or NA or NaN
  if (optim_args$method == "L-BFGS-B" || optim_args$method == "Brent") {
    big_finite_val <- 10 ^ 10
    neg_prof_loglik <- function(x) {
      pars <- numeric(p)
      pars[prof_pars] <- prof_vals
      pars[rest_pars] <- x
      check <- -do.call(object, list(pars, type = type))
      if (!is.finite(check)) {
        check <- big_finite_val
      }
      return(check)
    }
  }
  #
  for_optim <- c(list(par = init, fn = neg_prof_loglik), optim_args)
  #
  temp <- do.call(stats::optim, for_optim)
  #
  return(temp$value)
}
