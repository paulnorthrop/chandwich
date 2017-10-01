# Check length of init and fixed_pars are consistent

# ============================== compare_models  ==============================
#
#' Comparison of nested models
#'
#' adjusted likelihood ratio statistics (ALRS),
#' following Section 3.5 of
#' \href{http://dx.doi.org/10.1093/biomet/asm015}{Chandler and Bate (2007)}.
#'
#' @param larger An object of class \code{"chandwich"} returned by
#'   \code{adjust_loglik}.  The larger of the two models.
#' @param smaller An object of class \code{"chandwich"} returned by
#'   \code{adjust_loglik}.  The smaller of the two models.
#' @param approx A logical scalar.  If \code{approx = TRUE} then the
#'   approximation detailed by equations (18)-(20) of Chandler and Bates(2007)
#'   is used.  This option is only used if \code{smaller} is supplied.
#'   If \code{approx = FALSE} then the adjusted likelihood is
#'   maximised under the restriction imposed by \code{delta}.
#' @param fixed_pars A numeric vector.
#' @param fixed_at A numeric vector.
#' @param init A numeric vector of initial values for use in the search for
#'   the MLE under the smaller model.
#' @param ... Further arguments to be passed to \code{\link[stats]{optim}}.
#'   These may include \code{gr}, \code{method}, \code{lower}, \code{upper}
#'   or \code{control}.
#' @details add details
#' @return Test statistic, degrees of freedom (q) and p-value.
#' @references Chandler, R. E. and Bate, S. (2007). Inference for clustered
#'   data using the independence loglikelihood. \emph{Biometrika},
#'   \strong{94}(1), 167-183. \url{http://dx.doi.org/10.1093/biomet/asm015}
#' @seealso \code{\link{adjust_loglik}}: to adjust a user-supplied
#'   loglikelhood function.
#' @seealso \code{\link{summary.chandwich}} for maximum likelihood estimates
#'   and unadjusted and adjusted standard errors.
#' @seealso \code{\link{plot.chandwich}} for one- and two- dimensional plots
#'   of of adjusted loglikelihoods.
#' @examples
#' # GEV model, owtemps data ----------
#' # ... following Section 5.2 of Chandler and Bate (2007)
#'
#'   gev_loglik <- function(pars, data) {
#'     o_pars <- pars[c(1, 3, 5)] + pars[c(2, 4, 6)]
#'     w_pars <- pars[c(1, 3, 5)] - pars[c(2, 4, 6)]
#'     if (o_pars[2] <= 0 | w_pars[2] <= 0) return(-Inf)
#'     o_loglik <- revdbayes::dgev(data[, "Oxford"], o_pars[1], o_pars[2],
#'                                 o_pars[3], log = TRUE)
#'     w_loglik <- revdbayes::dgev(data[, "Worthing"], w_pars[1], w_pars[2],
#'                                 w_pars[3], log = TRUE)
#'     return(o_loglik + w_loglik)
#'   }
#'
#' # Initial estimates (method of moments for the Gumbel case)
#' sigma <- as.numeric(sqrt(6 * diag(stats::var(owtemps))) / pi)
#' mu <- as.numeric(colMeans(owtemps) - 0.57722 * sigma)
#' init <- c(mean(mu), -diff(mu) / 2, mean(sigma), -diff(sigma) / 2, 0, 0)
#' # Perform the log-likelihood adjustment
#' larger <- adjust_loglik(gev_loglik, data = owtemps, init = init,
#'           par_names = c("mu0", "mu1", "sigma0", "sigma1", "xi0", "xi1"))
#' # Rows 1, 3 and 4 of Table 2 of Chandler and Bate (2007)
#' round(attr(larger, "MLE"), 4)
#' round(attr(larger, "SE"), 4)
#' round(attr(larger, "adjSE"), 4)
#'
#' smaller <- adjust_loglik(gev_loglik, data = owtemps, init = init,
#'            par_names = c("mu0", "mu1", "sigma0", "sigma1", "xi0", "xi1"),
#'            fixed_pars = 6)
#'
#' compare_models(larger, smaller)
#' compare_models(larger, fixed_pars = 6)
#' compare_models(larger, smaller, approx = TRUE)
#' @export
compare_models <- function(larger, smaller = NULL, approx = FALSE,
                           fixed_pars = NULL, fixed_at = 0, init = NULL, ...) {
  #
  # Setup and checks -----------------------------------------------------------
  #
  # If smaller is supplied then .......
  if (!is.null(smaller)) {
    fixed_pars <- attr(smaller, "fixed_pars")
    fixed_at <- attr(smaller, "fixed_at")
    qq <- length(fixed_pars)
    p <- attr(larger, "p_current")
    pars <- numeric(p)
    s_mle <- attr(smaller, "MLE")
    l_mle <- attr(larger, "MLE")
    pars[fixed_pars] <- fixed_at
    free_pars <- (1:p)[-fixed_pars]
    pars[free_pars] <- s_mle
    max_loglik_smaller <- sum(do.call(larger, list(pars)))
    if (approx) {
      HA <- attr(larger, "HA")
      R <- solve(-HA)
      pjn <- solve(-R[fixed_pars, fixed_pars])
      num <- t(l_mle[fixed_pars] - fixed_at) %*% pjn %*%
        (l_mle[fixed_pars] - fixed_at)
      den <- t(l_mle - pars) %*% HA %*% (l_mle - pars)
      const <- num / den
      alrts <- 2 * const * (attr(larger, "max_loglik") - max_loglik_smaller)
      return(list(alrts = alrts, df = qq,
                  p_value = 1 - stats::pchisq(alrts, qq),
                  larger_mle = l_mle, smaller_mle = s_mle))
    } else {
      if (is.null(init)) {
        init <- attr(smaller, "MLE")
      }
    }
  }
  if (is.null(fixed_pars)) {
    stop("'fixed_pars' must be supplied")
  }
  # The number of parameters in the larger model
#  p <- length(attr(larger, "MLE"))
  p <- attr(larger, "p_current")
  # The number of parameters that are fixed
  qq <- length(fixed_pars)
  if (qq >= p) {
    stop("length(fixed_pars) must be smaller than attr(larger, ''MLE'')")
  }
  if (!(length(fixed_at) %in% c(1, qq))) {
    stop("the lengths of 'fixed_pars' and 'fixed_at' are not compatible")
  }
  fixed_at <- rep_len(fixed_at, qq)
  free_pars <- (1:p)[-fixed_pars]
  #
  # Extract arguments to be passed to optim()
  optim_args <- list(...)
  # Use "BFGS", unless the user has chosen the method or if p = 1 and they
  # have (inappropriately) chosen "Nelder-Mead" when p = 1
  if (is.null(optim_args$method)) {
    optim_args$method <- "BFGS"
  } else if (p == 1 & optim_args$method == "Nelder-Mead") {
    optim_args$method <- "BFGS"
  }
  # If approx = TRUE then we maximise the independence loglikelihood
  # under the constraint that certain parameter(s) are fixed.
  #
  # If approx = FALSE then we maximise the adjusted loglikelihood
  # under the constraint that certain parameter(s) are fixed.
  # Function to minimise to find restricted MLE of adjusted loglikelihood
  neg_adj_loglik <- function(x) {
    pars <- numeric(p)
    pars[fixed_pars] <- fixed_at
    pars[free_pars] <- x
    loglik_vals <- do.call(larger, list(pars))
    return(-sum(loglik_vals))
  }
  # L-BFGS-B and Brent don't like Inf or NA or NaN
  if (optim_args$method == "L-BFGS-B" || optim_args$method == "Brent") {
    big_finite_val <- 10 ^ 10
    neg_adj_loglik <- function(x) {
      pars <- numeric(p)
      pars[fixed_pars] <- fixed_at
      pars[free_pars] <- x
      loglik_vals <- do.call(larger, list(pars))
      check <- -sum(loglik_vals)
      if (!is.finite(check)) {
        check <- big_finite_val
      }
      return(check)
    }
  }
  # Initial estimate: the MLE with fixed_pars set at the values in fixed_at
  if (is.null(init)) {
    init <- attr(larger, "MLE")[-fixed_pars]
  }
  for_optim <- c(list(par = init, fn = neg_adj_loglik), optim_args)
  #
  temp <- do.call(stats::optim, for_optim)
  alrts <- 2 * (attr(larger, "max_loglik") + temp$value)
  return(list(alrts = alrts, df = qq, p_value = 1 - stats::pchisq(alrts, qq),
              larger_mle = attr(larger, "MLE"), smaller_mle = temp$par))
}
