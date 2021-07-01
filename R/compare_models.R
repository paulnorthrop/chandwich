# ============================== compare_models ===============================

#' Comparison of nested models
#'
#' Compares nested models using the adjusted likelihood ratio test statistic
#' (ALRTS) described in Section 3.5 of Chandler and Bate (2007). The nesting
#' must result from the simple constraint that a subset of the parameters of
#' the larger model is held fixed.
#'
#' @param larger An object of class \code{"chandwich"} returned by
#'   \code{\link{adjust_loglik}}.  The larger of the two models.
#' @param smaller An object of class \code{"chandwich"} returned by
#'   \code{\link{adjust_loglik}}.  The smaller of the two models.
#'
#'   If \code{smaller} is supplied then the arguments \code{fixed_pars} and
#'   \code{fixed_at} described below are ignored.
#' @param approx A logical scalar.  If \code{approx = TRUE} then the
#'   approximation detailed by equations (18)-(20) of Chandler and Bate (2007)
#'   is used.  This option is available only if \code{smaller} is supplied.
#'   If \code{smaller} is not supplied then \code{approx = TRUE} is used,
#'   with no warning.
#'
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
#'   For full details see Section 3.5 of Chandler and Bate (2007).
#'   If \code{approx = FALSE} then the a likelihood ratio test of the null
#'   hypothesis that the smaller model is a valid simplification of the larger
#'   model is carried out directly using equation (17) of Chandler and Bate
#'   (2007) based on the adjusted loglikelihood under the larger model,
#'   returned by \code{adjust_loglik}.  This adjusted loglikelihood is
#'   maximised subject to the constraint that a subset of the parameters
#'   in the larger model are fixed.  If \code{smaller} is supplied
#'   then this maximisation can be avoided using an approximation
#'   detailed by equations (18)-(20) of Chandler and Bate (2007), which uses
#'   the MLE under the smaller model.  The same null distribution (chi-squared
#'   with degrees of freedom equal to the number of parameters that are fixed)
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
#'   \strong{94}(1), 167-183. \doi{10.1093/biomet/asm015}
#' @seealso \code{\link{adjust_loglik}} to adjust a user-supplied
#'   loglikelihood function.
#' @seealso \code{\link{conf_intervals}} for confidence intervals for
#'   individual parameters.
#' @seealso \code{\link{conf_region}} for a confidence region for
#'   pairs of parameters.
#' @seealso \code{\link{print.compmod}}.
#' @examples
#' # -------------------------- GEV model, owtemps data -----------------------
#' # ------------ following Section 5.2 of Chandler and Bate (2007) -----------
#'
#' gev_loglik <- function(pars, data) {
#'   o_pars <- pars[c(1, 3, 5)] + pars[c(2, 4, 6)]
#'   w_pars <- pars[c(1, 3, 5)] - pars[c(2, 4, 6)]
#'   if (isTRUE(o_pars[2] <= 0 | w_pars[2] <= 0)) return(-Inf)
#'   o_data <- data[, "Oxford"]
#'   w_data <- data[, "Worthing"]
#'   check <- 1 + o_pars[3] * (o_data - o_pars[1]) / o_pars[2]
#'   if (isTRUE(any(check <= 0))) return(-Inf)
#'   check <- 1 + w_pars[3] * (w_data - w_pars[1]) / w_pars[2]
#'   if (isTRUE(any(check <= 0))) return(-Inf)
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
#' compare_models(large, fixed_pars = "xi[1]")
#' compare_models(large, medium)
#' # Test xi1 = 0, using approximation
#' compare_models(large, medium, approx = TRUE)
#'
#' # Horizontal adjustments
#' compare_models(large, medium, type = "cholesky")$p_value
#' compare_models(large, medium, type = "spectral")$p_value
#' # No adjustment (independence loglikelihood)
#' compare_models(large, medium, type = "none")$p_value
#'
#' # Test sigma1 = 0 for model with xi1 = 0
#' compare_models(medium, small)
#' # Test sigma1 = xi1 = 0
#' compare_models(large, small)
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
  p_l <- attr(larger, "p_current")
  # The number of parameters in the full model
  p_full <- attr(larger, "p_full")
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
    nest_1 <- p_l > p_s
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
    qq <- p_l - p_s
    s_mle <- attr(smaller, "MLE")
    l_mle <- attr(larger, "MLE")
    # To ensure that the parameter vectors passed to the models have the
    # correct parameter values in the correct places first set up a
    # parameter vector with length equal to the number of parameters in the
    # full model, allocate values, and then prune if necessary
    pars <- numeric(p_full)
    pars[fixed_pars] <- fixed_at
    free_pars <- (1:p_full)[-fixed_pars]
    pars[free_pars] <- s_mle
    if (!is.null(attr(larger, "fixed_pars"))) {
      pars <- pars[-attr(larger, "fixed_pars")]
    }
    names(pars) <- attr(larger, "par_names")
    if (approx) {
      max_loglik_smaller <- do.call(larger, list(pars, type = type))
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
  } else {
    # If smaller is not supplied and fixed_pars is numeric then infer the names
    # of the fixed parameters, if these are available in larger
    if (!is.null(attr(larger, "full_par_names"))) {
      names(fixed_pars) <- attr(larger, "full_par_names")[fixed_pars]
    }
  }
  if (is.null(fixed_pars)) {
    stop("'fixed_pars' must be supplied")
  } else {
    # The number of parameters that are fixed
    len_fixed_pars <- length(fixed_pars)
    fixed_at <- rep_len(fixed_at, len_fixed_pars)
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
  if (len_fixed_pars >= p_l) {
    stop("length(fixed_pars) must be smaller than attr(larger, ''MLE'')")
  }
  if (!(length(fixed_at) %in% c(1, len_fixed_pars))) {
    stop("the lengths of 'fixed_pars' and 'fixed_at' are not compatible")
  }
  free_pars <- (1:p_full)[-fixed_pars]
  p_s <- length(free_pars)
  qq <- p_l - p_s
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
    pars <- numeric(p_full)
    pars[fixed_pars] <- fixed_at
    pars[free_pars] <- x
    if (!is.null(attr(larger, "fixed_pars"))) {
      pars <- pars[-attr(larger, "fixed_pars")]
    }
    names(pars) <- attr(larger, "par_names")
    loglik_vals <- do.call(larger, list(pars, type = type))
    return(-loglik_vals)
  }
  # L-BFGS-B and Brent don't like Inf or NA or NaN
  if (optim_args$method == "L-BFGS-B" || optim_args$method == "Brent") {
    big_finite_val <- 10 ^ 10
    neg_adj_loglik <- function(x) {
      pars <- numeric(p_full)
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

# ============================= anova.chandwich ===============================

#' Comparison of nested models
#'
#' \code{anova} method for objects of class \code{"chandwich"}.
#' Compares two or more nested models using the adjusted likelihood ratio
#' test statistic (ALRTS) described in Section 3.5 of Chandler and Bate (2007).
#' The nesting must result from the simple constraint that a subset of the
#' parameters of the larger model is held fixed.
#'
#' @param object An object of class \code{"chandwich"}, returned by
#'   \code{\link{adjust_loglik}}.
#' @param object2 An object of class \code{"chandwich"}, returned by
#'   \code{\link{adjust_loglik}}.
#' @param ... Further objects of class \code{"chandwich"} and/or arguments
#'   to be passed to \code{\link{compare_models}}.  The name of any object
#'   of class \code{"chandwich"} passed via ... must not match any argument of
#'   \code{\link{compare_models}} or any argument of
#'   \code{\link[stats]{optim}}.
#' @details For details the adjusted likelihood ratio test see
#' \code{\link{compare_models}} and Chandler and Bate (2007).
#'
#'   The objects of class \code{"chandwich"} need not be provided in nested
#'   order: they will be ordered inside \code{anova.chandwich} based on the
#'   values of \code{attr(., "p_current")}.
#' @return An object of class \code{"anova"} inheriting from class
#'  \code{"data.frame"}, with four columns:
#'     \item{Model.Df}{The number of parameters in the model}
#'     \item{Df}{The decrease in the number of parameter compared the model
#'       in the previous row}
#'     \item{ALRTS}{The adjusted likelihood ratio test statistic}
#'     \item{Pr(>ALRTS)}{The p-value associated with the test that the
#'       model is a valid simplification of the model in the previous row.}
#'  The row names are the names of the model objects.
#' @seealso \code{\link{compare_models}} for an adjusted likelihood ratio test
#'   of two models.
#' @seealso \code{\link{adjust_loglik}} to adjust a user-supplied
#'   loglikelihood function.
#' @seealso \code{\link{conf_intervals}} for confidence intervals for
#'   individual parameters.
#' @seealso \code{\link{conf_region}} for a confidence region for
#'   pairs of parameters.
#' @references Chandler, R. E. and Bate, S. (2007). Inference for clustered
#'   data using the independence loglikelihood. \emph{Biometrika},
#'   \strong{94}(1), 167-183. \doi{10.1093/biomet/asm015}
#' @examples
#' # -------------------------- GEV model, owtemps data -----------------------
#' # ------------ following Section 5.2 of Chandler and Bate (2007) -----------
#'
#' gev_loglik <- function(pars, data) {
#'   o_pars <- pars[c(1, 3, 5)] + pars[c(2, 4, 6)]
#'   w_pars <- pars[c(1, 3, 5)] - pars[c(2, 4, 6)]
#'   if (isTRUE(o_pars[2] <= 0 | w_pars[2] <= 0)) return(-Inf)
#'   o_data <- data[, "Oxford"]
#'   w_data <- data[, "Worthing"]
#'   check <- 1 + o_pars[3] * (o_data - o_pars[1]) / o_pars[2]
#'   if (isTRUE(any(check <= 0))) return(-Inf)
#'   check <- 1 + w_pars[3] * (w_data - w_pars[1]) / w_pars[2]
#'   if (isTRUE(any(check <= 0))) return(-Inf)
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
#' tiny <- adjust_loglik(larger = small,
#'                       fixed_pars = c("mu[1]", "sigma[1]", "xi[1]"))
#'
#' anova(large, medium, small, tiny)
#' @export
anova.chandwich <- function (object, object2, ...) {
  if (missing(object)) {
    stop("model one must be supplied, using object")
  }
  if (missing(object2)) {
    stop("model two must be supplied, using object2")
  }
  # Extract the names of object and object2
  model1 <- deparse(substitute(object))
  model2 <- deparse(substitute(object2))
  # Extract arguments supplied in ... and determine which are named
  dotargs <- list(...)
  named <- if (is.null(names(dotargs)))
    rep_len(FALSE, length(dotargs))
  else (names(dotargs) != "")
  which_named <- which(named)
  which_not_named <- which(!named)
  # Named objects are intended for compare_models()
  for_compare_models <- dotargs[named]
  # Create list of model objects:  unnamed arguments may be model objects
  model_list <- c(list(object, object2), dotargs[!named])
  # Check for objects that do not have class "chandwich"
  is_chand <- vapply(model_list, function(x) inherits(x, "chandwich"), NA)
  if (any(!is_chand)) {
    stop("The following are not 'chandwich' objects: ",
         paste(names(model_list)[!is_chand], collapse = ", "))
  }
  extra_names <- as.list(substitute(list(...)))[-1][which_not_named]
  extra_names <- sapply(extra_names, function(x) deparse(x))
  model_names <- c(model1, model2, extra_names)
  # Check for duplicate names
  if (anyDuplicated(model_names)) {
    stop("A model name has been supplied more than once")
  }
  # Order the models in order of the number of parameters
  n_pars <- vapply(model_list, function(x) attr(x, "p_current"), 0)
  # Check for models with the same number of parameters
  if (anyDuplicated(n_pars)) {
    stop("At least two models have the same number of parameters")
  }
  m_order <- order(n_pars, decreasing = TRUE)
  model_list <- model_list[m_order]
  n_pars <- n_pars[m_order]
  n_models <- length(model_list)
  # Do the testing
  alrts <- p_value <- numeric(n_models - 1)
  for (i in 2:n_models) {
    larger <- model_list[[i - 1]]
    smaller <- model_list[[i]]
    res <- do.call(compare_models, c(list(larger = larger, smaller = smaller),
                                     for_compare_models))
    alrts[i - 1] <- res$alrts
    p_value[i - 1] <- res$p_value
  }
  df <- -diff(n_pars)
  my_table <- data.frame(n_pars, c(NA, df), c(NA, alrts), c(NA, p_value))
  dimnames(my_table) <- list(model_names,
                             c("Model.Df", "Df", "ALRTS", "Pr(>ALRTS)"))
  structure(my_table, heading = c("Analysis of (Adjusted) Deviance Table\n"),
            class = c("anova", "data.frame"))
}
