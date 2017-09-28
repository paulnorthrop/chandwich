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
#' @param fixed_pars A numeric vector.
#' @param fixed_at A numeric vector.
#' @param init A numeric vector of initial values for use in the search for
#'   the MLE under the smaller model.
#' @param approx A logical scalar.  If \code{approx = TRUE} then the
#'   approximation detailed by equations (18)-(20) of Chandler and Bates(2007)
#'   is used.  If \code{approx = FALSE} then the adjusted likelihood is
#'   maximised under the restriction imposed by \code{delta}.
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
#' if (requireNamespace("revdbayes", quietly = TRUE)) {
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
#' }
#' # Initial estimates (method of moments for the Gumbel case)
#' sigma <- as.numeric(sqrt(6 * diag(stats::var(owtemps))) / pi)
#' mu <- as.numeric(colMeans(owtemps) - 0.57722 * sigma)
#' init <- c(mean(mu), -diff(mu) / 2, mean(sigma), -diff(sigma) / 2, 0, 0)
#' # Perform the log-likelihood adjustment
#' ow_res <- adjust_loglik(gev_loglik, data = owtemps, init = init,
#'           par_names = c("mu0", "mu1", "sigma0", "sigma1", "xi0", "xi1"))
#' # Rows 1, 3 and 4 of Table 2 of Chandler and Bate (2007)
#' round(attr(ow_res, "MLE"), 4)
#' round(attr(ow_res, "SE"), 4)
#' round(attr(ow_res, "adjSE"), 4)
#'
#' compare_models(ow_res, fixed_pars = 6)
#' @export
compare_models <- function(larger, smaller = NULL, fixed_pars = NULL,
                           fixed_at = 0, init = NULL, approx = TRUE, ...) {
  #
  # Setup and checks -----------------------------------------------------------
  #
  if (is.null(fixed_pars)) {
    stop("'fixed_pars' must be supplied")
  }
  # Indices and/or numbers ...................................
  qq <- length(fixed_pars)
  if (!(length(fixed_at) %in% c(1, qq))) {
    stop("the lengths of 'fixed_pars' and 'fixed_at' are not compatible")
  }
  fixed_at <- rep_len(fixed_at, qq)
  other_pars <- (1:attr(ow_res, "p"))[-fixed_pars]
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
  # If quick = FALSE then we need to maximise the adjusted loglikelihood
  # under the constraint that certain parameter(s) are fixed.
  #
  # Function to minimise to find restricted MLE of adjusted loglikelihood
  neg_adj_loglik <- function(x) {
    pars <- numeric(qq)
    pars[fixed_pars] <- fixed_at
    pars[other_pars] <- x
    loglik_vals <- do.call(larger, list(pars))
    return(-sum(loglik_vals))
  }
  # L-BFGS-B and Brent don't like Inf or NA or NaN
  if (optim_args$method == "L-BFGS-B" || optim_args$method == "Brent") {
    big_finite_val <- 10 ^ 10
    neg_adj_loglik <- function(x) {
      pars <- numeric(qq)
      pars[fixed_pars] <- fixed_at
      pars[other_pars] <- x
      loglik_vals <- do.call(larger, list(x))
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
  temp <- do.call(optim, for_optim)
  alrts <- 2 * (attr(larger, "max_loglik") + temp$value)
  return(list(alrts = alrts, df = qq, p_value = 1 - pchisq(alrts, qq),
              larger_mle = attr(larger, "MLE"), smaller_mle = temp$par))
}
