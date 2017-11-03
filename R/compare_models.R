# ============================== compare_models ===============================

#' Comparison of nested models
#'
#' Compares nested models using the adjusted likelihood ratio statistic (ALRS)
#' described in Section 3.5 of
#' \href{http://dx.doi.org/10.1093/biomet/asm015}{Chandler and Bate (2007)}.
#' The nesting must result from the simple constraint that a subset of the
#' parameters of the larger model is held fixed.
#'
#' @param larger An object of class \code{"chandwich"} returned by
#'   \code{adjust_loglik}.  The larger of the two models.
#' @param smaller An object of class \code{"chandwich"} returned by
#'   \code{adjust_loglik}.  The smaller of the two models.
#'
#'   If \code{smaller} is supplied then the arguments \code{fixed_pars} and
#'   \code{fixed_at} are ignored.
#' @param approx A logical scalar.  If \code{approx = TRUE} then the
#'   approximation detailed by equations (18)-(20) of Chandler and Bate (2007)
#'   is used.  This option is available only if \code{smaller} is supplied.
#'   The approximation doesn't make sense if \code{type = "none"}.  If
#'   \code{type = "none"} and \code{approx = TRUE} then \code{approx} is
#'   set to \code{FALSE} with no warning.
#'
#'   If \code{approx = FALSE} then the adjusted likelihood is
#'   maximised under the restriction imposed by \code{delta}, that is,
#'   equation (17) of Chandler and Bate (2007) is used.
#' @param type A character scalar.  The argument \code{type} to the function
#'   returned by \code{\link{adjust_loglik}}, that is, the type of adjustment
#'   made to the independence loglikelihood function.
#' @param fixed_pars A numeric vector.  Indices of the components of the
#'   \strong{full} parameter vector that are restricted to be equal to the
#'   value(s) in \code{fixed_at}.
#' @param fixed_at A numeric vector of length 1 or \code{length(fixed_pars)}.
#'   If \code{length(fixed_at) = 1} then the components \code{fixed_pars}
#'   of the parameter vector are all fixed at \code{fixed_at}.
#'   If \code{length(fixed_at) = length(fixed_pars)} then the component
#'   \code{fixed_pars[i]} is fixed at \code{fixed_at[i]} for each \code{i}.
#' @param init (Only relevant if \code{approx = FALSE}).
#'   A numeric vector of initial values for use in the search for
#'   the MLE under the smaller model.  Must have length equal to the number
#'   of parameters in the smaller of the two models being compared.  If
#'   \code{init} is not supplied, or it is of the wrong length, then
#'   \code{attr(smaller, "MLE")} is used if \code{smaller} is supplied and
#'   \code{attr(larger, "MLE")[-fixed_pars]} otherwise.
#' @param ... Further arguments to be passed to \code{\link[stats]{optim}}.
#'   These may include \code{gr}, \code{method}, \code{lower}, \code{upper}
#'   or \code{control}.
#' @details The smaller of the two models is specified either by supplying
#'   \code{smaller} or \code{fixed_pars}.  If both are supplied then
#'   \code{smaller} takes precedence.
#'
#'   For full details see Section 3.5 of
#'   \href{http://dx.doi.org/10.1093/biomet/asm015}{Chandler and Bate (2007)}.
#'   if \code{approx = FALSE} then the a likelihood ratio test of the null
#'   hypothesis that the smaller model is a valid simplication of the larger
#'   model is carried out directly using equation (17) of Chandler and Bate
#'   (2007) based on the adjusted loglikelihood under the larger model,
#'   returned by \code{adjust_loglik}.  This adjusted loglikelihood is
#'   maximised subject to the constraint that a subset of the parameters
#'   in the larger model are fixed.  If \code{smaller} is supplied
#'   then this maximisation can be avoided using an approximation
#'   detailed by (18)-(20) of Chandler and Bate (2007), which uses the
#'   MLE under the smaller model.  The same null distribution (chi-squared with
#'   degrees of freedom equal to the number of parameters that are fixed)
#'   is used in both cases.
#' @return An object of class "compmod", a list with components
#'  \item{alrts}{the adjusted likelihood ratio test statistic.}
#'  \item{df}{under the null hypothesis that the smaller model is a valid
#'    simplification of the larger model the adjusted likelihood ratio
#'    test statistic has a chi-squared distribution with \code{df}
#'    degrees of freedom.}
#'  \item{p_value}{the p-value associated with the test of the null
#'    hypothesis.}
#'  \item{larger_mle}{the MLE of the parameters under the larger model.}
#'  \item{smaller_mle}{the MLE of the parameters under the smaller model.}
#'  \item{larger_fixed_pars, smaller_fixed_pars}{Numeric vectors of the
#'    indices of parameters fixed in the larger and smaller models,
#'    respectively.}
#'  \item{larger_fixed_at, smaller_fixed_at}{Numeric vectors of the
#'    values at which the parameters in \code{larger_fixed_pars} and
#'    \code{smaller_fixed_pars} are fixed.}
#'  \item{approx}{the argument \code{approx}.}
#' @references Chandler, R. E. and Bate, S. (2007). Inference for clustered
#'   data using the independence loglikelihood. \emph{Biometrika},
#'   \strong{94}(1), 167-183. \url{http://dx.doi.org/10.1093/biomet/asm015}
#' @seealso \code{\link{adjust_loglik}} to adjust a user-supplied
#'   loglikelhood function.
#' @seealso \code{\link{summary.chandwich}} for maximum likelihood estimates
#'   and unadjusted and adjusted standard errors.
#' @seealso \code{\link{plot.chandwich}} for plots of one-dimensional adjusted
#'   loglikelihoods.
#' @seealso \code{\link{conf_intervals}} for confidence intervals for
#'   individual parameters.
#' @seealso \code{\link{conf_region}} for a confidence region for
#'   pairs of parameters.
#' @examples
#' # -------------------------- GEV model, owtemps data -----------------------
#' # ------------ following Section 5.2 of Chandler and Bate (2007) -----------
#'
#' gev_loglik <- function(pars, data) {
#'   o_pars <- pars[c(1, 3, 5)] + pars[c(2, 4, 6)]
#'   w_pars <- pars[c(1, 3, 5)] - pars[c(2, 4, 6)]
#'   if (o_pars[2] <= 0 | w_pars[2] <= 0) return(-Inf)
#'   o_data <- data[, "Oxford"]
#'   w_data <- data[, "Worthing"]
#'   check <- 1 + o_pars[3] * (o_data - o_pars[1]) / o_pars[2]
#'   if (any(check <= 0)) return(-Inf)
#'   check <- 1 + w_pars[3] * (w_data - w_pars[1]) / w_pars[2]
#'   if (any(check <= 0)) return(-Inf)
#'   o_loglik <- log_gev(o_data, o_pars[1], o_pars[2], o_pars[3])
#'   w_loglik <- log_gev(w_data, w_pars[1], w_pars[2], w_pars[3])
#'   return(o_loglik + w_loglik)
#' }
#'
#' # Initial estimates (method of moments for the Gumbel case)
#' sigma <- as.numeric(sqrt(6 * diag(var(owtemps))) / pi)
#' mu <- as.numeric(colMeans(owtemps) - 0.57722 * sigma)
#' init <- c(mean(mu), -diff(mu) / 2, mean(sigma), -diff(sigma) / 2, 0, 0)
#'
#' # Log-likelihood adjustment of the full model
#' par_names <- c("mu[0]", "mu[1]", "sigma[0]", "sigma[1]", "xi[0]", "xi[1]")
#' large <- adjust_loglik(gev_loglik, data = owtemps, init = init,
#'          par_names = par_names)
#'
#' # Log-likelihood adjustment of some smaller models: xi[1] = 0 etc
#'
#' medium <- adjust_loglik(larger = large, fixed_pars = "xi[1]")
#' small <- adjust_loglik(larger = medium, fixed_pars = c("sigma[1]", "xi[1]"))
#'
#' # Tests
#'
#' # Test xi1 = 0 (2 equivalent ways), vertical adjustment
#' compare_models(large, fixed_pars = "xi[1]")$p_value
#' compare_models(large, medium)$p_value
#' # Test xi1 = 0, using approximation
#' compare_models(large, medium, approx = TRUE)$p_value
#'
#' # Horizontal adjustments
#' compare_models(large, medium, type = "cholesky")$p_value
#' compare_models(large, medium, type = "spectral")$p_value
#' # No adjustment (independence loglikelihood)
#' compare_models(large, medium, type = "none")$p_value
#'
#' # Test sigma1 = 0 for model with xi1 = 0
#' compare_models(medium, small)$p_value
#' # Test sigma1 = xi1 = 0
#' compare_models(large, small)$p_value
#'
#' # --------- Misspecified Poisson model for negative binomial data ----------
#'
#' # ... following Section 5.1 of the "Object-Oriented Computation of Sandwich
#' # Estimators" vignette of the sandwich package
#' # https://cran.r-project.org/web/packages/sandwich/vignettes/sandwich-OOP.pdf
#'
#' # Simulate data
#' set.seed(123)
#' x <- rnorm(250)
#' y <- rnbinom(250, mu = exp(1 + x), size = 1)
#' # Fit misspecified Poisson model
#' fm_pois <- glm(y ~ x + I(x^2), family = poisson)
#' summary(fm_pois)$coefficients
#'
#' # Contributions to the independence loglikelihood
#' pois_glm_loglik <- function(pars, y, x) {
#'   log_mu <- pars[1] + pars[2] * x + pars[3] * x ^ 2
#'   return(dpois(y, lambda = exp(log_mu), log = TRUE))
#' }
#' pars <- c("alpha", "beta", "gamma")
#' pois_quad <- adjust_loglik(pois_glm_loglik, y = y, x = x, par_names = pars)
#' summary(pois_quad)
#'
#' pois_lin <- adjust_loglik(larger = pois_quad, fixed_pars = "gamma")
#'
#' # Test the significance of the quadratic term
#' compare_models(pois_quad, pois_lin)$p_value
#' compare_models(pois_quad, pois_lin, approx = TRUE)$p_value
#' @export
compare_models <- function(larger, smaller = NULL, approx = FALSE,
                           type = c("vertical", "cholesky", "spectral",
                                    "none"), fixed_pars = NULL,
                           fixed_at = rep_len(0, length(fixed_pars)),
                           init = NULL, ...) {
  type <- match.arg(type)
  if (type == "none" & approx) {
    approx <- FALSE
  }
  #
  # Setup and checks -----------------------------------------------------------
  #
  # The number of parameters in the larger model
  p <- attr(larger, "p_current")
  # Extract the fixed parameters (if any) from the larger model
  l_fixed_pars <- attr(larger, "fixed_pars")
  l_fixed_at <- attr(larger, "fixed_at")
  # If smaller is supplied then .......
  if (!is.null(smaller)) {
    # Check that smaller is nested within larger
    # (1) larger and smaller must be derived from the same full model
    if (attr(larger, "name") != attr(smaller, "name")) {
      stop("larger and smaller are not derived from the same model")
    }
    # (2) larger must have more parameters than smaller
    # (3) all values in attr(larger, "fixed_pars") must also appear in
    #     attr(smaller, "fixed_pars")
    # (4) Any parameters that appear in both attr(larger, "fixed_pars") and
    #     attr(smaller, "fixed_pars") must have the same corresponding values
    #     in attr(larger, "fixed_at") and attr(smaller, "fixed_at")
    p_s <- attr(smaller, "p_current")
    nest_1 <- p > p_s
    fixed_pars <- attr(smaller, "fixed_pars")
    fixed_at <- attr(smaller, "fixed_at")
    s_fixed_pars <- fixed_pars
    s_fixed_at <- fixed_at
    nest_2 <- all(l_fixed_pars %in% s_fixed_pars)
    if (!nest_1 | !nest_2) {
      stop("smaller is not nested in larger")
    }
    l_in_s <- which(l_fixed_pars %in% s_fixed_pars)
    s_in_l <- which(s_fixed_pars %in% l_fixed_pars)
    if (any(l_fixed_at[l_in_s] != s_fixed_at[s_in_l])) {
      stop("smaller is not nested in larger: ",
           "parameter(s) fixed at different values")
    }
    qq <- length(fixed_pars)
    p <- attr(larger, "p_current")
    s_mle <- attr(smaller, "MLE")
    l_mle <- attr(larger, "MLE")
    pars <- numeric(p)
    pars[fixed_pars] <- fixed_at
    free_pars <- (1:p)[-fixed_pars]
    pars[free_pars] <- s_mle
    if (!is.null(attr(larger, "fixed_pars"))) {
      pars <- pars[-attr(larger, "fixed_pars")]
    }
    max_loglik_smaller <- do.call(larger, list(pars, type = type))
    if (approx) {
      HA <- attr(larger, "HA")
      R <- solve(-HA)
      # We want only the parameters that are fixed in smaller but not in larger
      s_not_in_l <- which(!(s_fixed_pars %in% l_fixed_pars))
      new_fixed_pars <- s_fixed_pars[s_not_in_l]
      new_fixed_at <- fixed_at[s_not_in_l]
      R_rest <- solve(-R[new_fixed_pars, new_fixed_pars])
      num <- t(l_mle[new_fixed_pars] - new_fixed_at) %*% R_rest %*%
        (l_mle[new_fixed_pars] - new_fixed_at)
      den <- t(l_mle - pars) %*% HA %*% (l_mle - pars)
      const <- num / den
      alrts <- 2 * const * (attr(larger, "max_loglik") - max_loglik_smaller)
      alrts <- as.numeric(alrts)
      comp_list <- list(alrts = alrts, df = qq,
                        p_value = 1 - stats::pchisq(alrts, qq),
                        larger_mle = l_mle, smaller_mle = s_mle,
                        name = attr(larger, "name"),
                        larger_fixed_pars = l_fixed_pars,
                        larger_fixed_at = l_fixed_at,
                        smaller_fixed_pars = s_fixed_pars,
                        smaller_fixed_at = s_fixed_at, approx = approx)
      class(comp_list) <- "compmod"
      return(comp_list)
    } else {
      if (is.null(init) || length(init) != p_s) {
        init <- attr(smaller, "MLE")
      }
    }
  }
  if (is.null(fixed_pars)) {
    stop("'fixed_pars' must be supplied")
  } else {
    fixed_at <- rep_len(fixed_at, length(fixed_pars))
    # If fixed_pars is a character vector then
    # (a) check that full_par_names is not NULL
    # (b) check that fixed_pars is a subset of full_par_names
    # (c) determine the parameter indices of the components of fixed_pars
    if (is.character(fixed_pars)) {
      if (is.null(attr(larger, "full_par_names"))) {
        stop("fixed_pars can be character only if par_names is supplied")
      }
      full_par_names <- attr(larger, "full_par_names")
      if (!all(fixed_pars %in% full_par_names)) {
        stop("fixed_pars is not a subset of ", deparse(full_par_names))
      }
      temp <- fixed_pars
      fixed_pars <- which(full_par_names %in% fixed_pars)
      names(fixed_pars) <- temp
    }
  }
  s_fixed_pars <- fixed_pars
  s_fixed_at <- fixed_at
  # Check that the implied smaller model is nested within larger
  # (1) all values in l_fixed_pars must also appear in s_fixed_pars
  if (!all(l_fixed_pars %in% s_fixed_pars)) {
    stop("smaller is not nested in larger")
  }
  # (2) Any parameters that appear in both attr(larger, "fixed_pars") and
  #     attr(smaller, "fixed_pars") must have the same corresponding values
  #     in attr(larger, "fixed_at") and attr(smaller, "fixed_at")
  l_in_s <- which(l_fixed_pars %in% s_fixed_pars)
  s_in_l <- which(s_fixed_pars %in% l_fixed_pars)
  if (any(l_fixed_at[l_in_s] != s_fixed_at[s_in_l])) {
    stop("smaller is not nested in larger: ",
         "parameter(s) fixed at different values")
  }
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
  p_s <- length(free_pars)
  #
  # Extract arguments to be passed to optim()
  optim_args <- list(...)
  # Use "BFGS", unless the user has chosen the method or if p_s = 1 and they
  # have (inappropriately) chosen "Nelder-Mead" when p_s = 1
  if (is.null(optim_args$method)) {
    optim_args$method <- "BFGS"
  } else if (p_s == 1 & optim_args$method == "Nelder-Mead") {
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
    if (!is.null(attr(larger, "fixed_pars"))) {
      pars <- pars[-attr(larger, "fixed_pars")]
    }
    loglik_vals <- do.call(larger, list(pars, type = type))
    return(-loglik_vals)
  }
  # L-BFGS-B and Brent don't like Inf or NA or NaN
  if (optim_args$method == "L-BFGS-B" || optim_args$method == "Brent") {
    big_finite_val <- 10 ^ 10
    neg_adj_loglik <- function(x) {
      pars <- numeric(p)
      pars[fixed_pars] <- fixed_at
      pars[free_pars] <- x
      if (!is.null(attr(larger, "fixed_pars"))) {
        pars <- pars[-attr(larger, "fixed_pars")]
      }
      loglik_vals <- do.call(larger, list(pars, type = type))
      check <- -loglik_vals
      if (!is.finite(check)) {
        check <- big_finite_val
      }
      return(check)
    }
  }
  # Initial estimate: the MLE with fixed_pars set at the values in fixed_at
  if (is.null(init) || length(init) != p_s) {
    init <- attr(larger, "MLE")[-fixed_pars]
  }
  for_optim <- c(list(par = init, fn = neg_adj_loglik), optim_args)
  #
  temp <- do.call(stats::optim, for_optim)
  alrts <- 2 * (attr(larger, "max_loglik") + temp$value)
  comp_list <- list(alrts = alrts, df = qq,
                    p_value = 1 - stats::pchisq(alrts, qq),
                    larger_mle = attr(larger, "MLE"), smaller_mle = temp$par,
                    name = attr(larger, "name"),
                    larger_fixed_pars = l_fixed_pars,
                    larger_fixed_at = l_fixed_at,
                    smaller_fixed_pars = s_fixed_pars,
                    smaller_fixed_at = s_fixed_at, approx = approx)
  class(comp_list) <- "compmod"
  return(comp_list)
}
