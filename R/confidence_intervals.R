# ============================== conf_region ===============================

#' Two-dimensional confidence regions
#'
#' Calculates the (profile, if necessary) loglikelihood for a pair of
#' parameters from which confidence regions can be plotted using
#' \code{\link{plot.confreg}}.
#'
#' @param object An object of class \code{"chandwich"} returned by
#'   \code{adjust_loglik}.
#' @param which_pars A vector of length 2 specifying the 2 (unfixed)
#'   parameters for which confidence region is required.
#'   Can be either a numeric vector, specifying indices of the components
#'   of the \strong{full} parameter vector, or a character vector of
#'   parameter names, which must be a subset of those supplied in
#'   \code{par_names} in the call to \code{\link{adjust_loglik}} that
#'   produced \code{object}.
#'
#'   \code{which_pars} must not have any parameters in common with
#'   \code{attr(object, "fixed_pars")}.  \code{which_pars} must not contain
#'   all of the unfixed parameters, i.e. there is no point in profiling over
#'   all the unfixed parameters.
#'
#'   If \code{which_pars} is not supplied but the current model has exactly
#'   two free parameters, i.e. \code{attr(object, "p_current") = 2} then
#'   \code{which_pars} is set to \code{attr(object, "free_pars") = 2}.
#' @param range1,range2 Numeric vectors of length 2.  Respective ranges (of
#'   the form \code{lower, upper}) of values of \code{which_pars[1]} and
#'   \code{which_pars[2]} over which to profile.
#'   Missing values in \code{range1} and/or \code{range2} are
#'   filled in using \code{conf} and \code{mult}.  See below for details.
#' @param conf A numeric scalar in (0, 100).  The highest confidence level
#'   of interest. This is only relevant if \code{lower} and \code{upper} are
#'   not supplied.  In that event \code{conf} is used, in combination with
#'   \code{mult}, to try to set up the grid of parameter values to include
#'   the largest confidence region of interest.
#' @param mult A numeric vector of length 1 or the same length as
#'   \code{which_pars}.
#'   The search for the profile loglikelihood-based confidence limits is
#'   conducted over the corresponding symmetric confidence intervals, extended
#'   by a factor of the corresponding component of \code{mult}.
#' @param num A numeric vector of length 1 or 2.  The numbers of values at which
#'   to evaluate the profile loglikelihood either side of the MLE.
#'   \code{num[i]} relates to \code{which_pars[i]}.  If \code{num} has length
#'   1 then \code{num} is replicated to have length 2.
#' @param type A character scalar.  The argument \code{type} to the function
#'   returned by \code{\link{adjust_loglik}}, that is, the type of adjustment
#'   made to the independence loglikelihood function.
#' @param ... Further arguments to be passed to \code{\link[stats]{optim}}.
#'   These may include \code{gr}, \code{method}, \code{lower}, \code{upper}
#'   or \code{control}.  Any arguments that are not appropriate for
#'   \code{\link[stats]{optim}}, i.e. not in
#'   \code{methods::formalArgs(stats::optim)},
#'   will be removed without warning.
#' @return An object of class "confreg", a list with components
#'     \item{grid1, grid2}{Numeric vectors.   Respective values of
#'       \code{which_pars[1]} and \code{which_pars[2]} in the grid over which
#'       the (profile) loglikelihood is evaluated. }
#'     \item{max_loglik}{A numeric scalar.  The value value of
#'       the loglikelihood at its maxium.}
#'     \item{prof_loglik}{An 2 \code{num} + 1 by 2 \code{num} + 1
#'       numeric matrix containing the values of the (profile) loglikelihood.}
#'     \item{type}{A character scalar. The input \code{type}.}
#'     \item{which_pars}{A numeric or character vector.  The input
#'       \code{which_pars}.  If the \code{which_pars} was numeric then
#'       it is supplemented by the parameter names, if these are available
#'       in \code{object}.}
#'     \item{name}{A character scalar. The name of the model,
#'       stored in \code{attr(object, "name")}.}
#' @seealso \code{\link{adjust_loglik}} to adjust a user-supplied
#'   loglikelihood function.
#' @seealso \code{\link{conf_intervals}} for confidence intervals for
#'   individual parameters.
#' @seealso \code{\link{compare_models}} to compare nested models using an
#'   (adjusted) likelihood ratio test.
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
#'                        par_names = par_names)
#'
#' \dontrun{
#' # Plots like those in Figure 4 of Chandler and Bate (2007)
#' # (a)
#' which_pars <- c("mu[0]", "mu[1]")
#' reg_1 <- conf_region(large, which_pars = which_pars)
#' reg_none_1 <- conf_region(large, which_pars = which_pars, type = "none")
#' plot(reg_1, reg_none_1)
#' # (b)
#' which_pars <- c("sigma[0]", "sigma[1]")
#' reg_2 <- conf_region(large, which_pars = which_pars)
#' reg_none_2 <- conf_region(large, which_pars = which_pars, type = "none")
#' plot(reg_2, reg_none_2)
#' # (c)
#' # Note: the naive and bivariate model contours are the reversed in the paper
#' which_pars <- c("sigma[0]", "xi[0]")
#' reg_3 <- conf_region(large, which_pars = which_pars)
#' reg_none_3 <- conf_region(large, which_pars = which_pars, type = "none")
#' plot(reg_3, reg_none_3)
#' }
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
#' # Linear model (gamma fixed at 0)
#' pois_lin <- adjust_loglik(pois_glm_loglik, y = y, x = x, par_names = pars,
#'                           fixed_pars = "gamma")
#' pois_vertical <- conf_region(pois_lin)
#' pois_none <- conf_region(pois_lin, type = "none")
#' plot(pois_none, pois_vertical, conf = c(50, 75, 95, 99), col = 2:1, lwd = 2,
#'      lty = 1)
#' @export
conf_region <- function(object, which_pars = NULL, range1 = c(NA, NA),
                        range2 = c(NA, NA), conf = 95, mult = 2, num = c(10, 10),
                        type = c("vertical", "cholesky", "spectral", "none"),
                        ...) {
  type <- match.arg(type)
  # Check that arguments supplied in ... can be passed to stats::optim()
  optim_args <- list(...)
  optim_names <- names(optim_args)
  ok <- optim_names %in% methods::formalArgs(stats::optim)
  optim_args <- optim_args[ok]
  # Adjust conf to make it more applicable to the marginal intervals used to
  # set the default grid on which the profile loglikelihood is calculated
  conf_for_search <- sqrt(conf) * 10
  num <- rep_len(num, 2)
  # If which_pars has not been supplied and there are 2 free parameters in
  # the current model then set which_pars to these 2 parameters
  if (is.null(which_pars) & attr(object, "p_current") == 2) {
    which_pars <- attr(object, "free_pars")
  }
  n_which_pars <- length(which_pars)
  if (n_which_pars != 2) {
    stop("which_pars must be a vector of length 2")
  }
  if (length(range1) != 2) {
    stop("range1 must be a vector of length 2")
  }
  if (length(range2) != 2) {
    stop("range2 must be a vector of length 2")
  }
  if (all(!is.na(range1))) {
    range1 <- sort(range1)
  }
  if (all(!is.na(range2))) {
    range2 <- sort(range2)
  }
  n_mult <- length(mult)
  if (n_mult != 1 & n_mult != n_which_pars) {
    stop("mult must have length 1 or the same length as which_pars")
  } else if (n_mult == 1) {
    mult <- rep(mult, n_which_pars)
  }
  # Fixed parameters, values at which they are fixed and parameter names
  fixed_pars <- attr(object, "fixed_pars")
  fixed_at <- attr(object, "fixed_at")
  full_par_names <- attr(object, "full_par_names")
  p <- attr(object, "p_full")
  # If which_par is a character
  if (is.character(which_pars)) {
    if (is.null(full_par_names)) {
      stop("which_pars can be character only if par_names is supplied")
    }
    if (!all(which_pars %in% full_par_names)) {
      stop("which_pars is not a subset of ", deparse(full_par_names))
    }
    temp <- which_pars
    which_pars <- which(full_par_names %in% which_pars)
    names(which_pars) <- temp
  } else {
    if (!is.null(full_par_names)) {
      names(which_pars) <- full_par_names[which_pars]
    }
  }
  if (any(which_pars %in% fixed_pars)) {
    stop("which_pars & attr(object,''fixed_pars'') have parameters in common")
  }
  # Marginal symmetric confidence intervals for each of the 2 parameters
  if (is.null(fixed_pars)) {
    free_pars <- 1:p
  } else {
    free_pars <- (1:p)[-fixed_pars]
  }
  res_mle <- attr(object, "res_MLE")
  res_SE <- res_adjSE <- numeric(p)
  res_SE[free_pars] <- attr(object, "SE")
  res_adjSE[free_pars] <- attr(object, "adjSE")
  z_val <- stats::qnorm(1 - (1 - conf_for_search / 100) / 2)
  which_mle <- res_mle[which_pars]
  #
  lower <- c(range1[1], range2[1])
  upper <- c(range1[2], range2[2])
  if (type == "none") {
    the_SEs <- res_SE[which_pars]
  } else {
    the_SEs <- res_adjSE[which_pars]
  }
  sym_lower <- which_mle - z_val * the_SEs
  sym_upper <- which_mle + z_val * the_SEs
  if (!is.finite(lower[1])) {
    lower[1] <- which_mle[1] - mult[1] * z_val * the_SEs[1]
  }
  if (!is.finite(upper[1])) {
    upper[1] <- which_mle[1] + mult[1] * z_val * the_SEs[1]
  }
  if (!is.finite(lower[2])) {
    lower[2] <- which_mle[2] - mult[2] * z_val * the_SEs[2]
  }
  if (!is.finite(upper[2])) {
    upper[2] <- which_mle[2] + mult[2] * z_val * the_SEs[2]
  }
  sym_CI <- cbind(sym_lower, sym_upper)
  colnames(sym_CI) <- c("lower", "upper")
  #
  # 2D profile loglikelihood
  #
  # Set up a grid of values of the parameters in which_pars
  temp1 <- seq(lower[1], which_mle[1], length.out = num[1] + 1)
  temp2 <- seq(which_mle[1], upper[1], length.out = num[1] + 1)[-1]
  grid1 <- c(temp1, temp2)
  temp1 <- seq(lower[2], which_mle[2], length.out = num[2] + 1)
  temp2 <- seq(which_mle[2], upper[2], length.out = num[2] + 1)[-1]
  grid2 <- c(temp1, temp2)
  leng1 <- length(grid1)
  leng2 <- length(grid2)
  # Matrix to store the profile loglikelihood values
  z <- matrix(nrow = leng1, ncol = leng2)
  # MLE of all free parameters, i.e. excluding fixed parameters and which_pars
  theta_start <- res_mle[-c(fixed_pars, which_pars)]
  #
  # To avoid illegal starting values, move from the MLE downwards and
  # then back up, at each stage using neighbouring estimates as
  # starting values.
  #
  # In which positions are the MLEs in the grid?
  mle_idx1 <- num[1] + 1
  mle_idx2 <- num[2] + 1
  eval_order1 <- c(rev(1:(mle_idx1 - 1)), mle_idx1:(2 * num[1] + 1))
  eval_order2 <- c(rev(1:(mle_idx2 - 1)), mle_idx2:(2 * num[2] + 1))
  start_array <- array(dim = c(length(theta_start), leng1, leng2))
  start_array[, mle_idx1, mle_idx2] <- theta_start
  for (i in eval_order1) {
    i_nbr <- max(i - 1, 1):(i + 1)
    i_nbr <- i_nbr[i_nbr <= leng1]
    for (j in eval_order2) {
      j_nbr <- max(j - 1, 1):(j + 1)
      j_nbr <- j_nbr[j_nbr <= leng2]
      theta <- rowMeans(start_array[, i_nbr, j_nbr], na.rm = TRUE)
      if (any(is.na(theta))) {
        theta <- theta_start
      }
      prof_vals <- c(grid1[i], grid2[j])
      prof_args <- c(list(object = object, prof_pars = which_pars,
                          prof_vals = prof_vals, init = theta, type = type),
                     optim_args)
      zz <- try(do.call(profile_loglik, prof_args), silent = TRUE)
      if (class(zz) == "try-error") {
        z[i, j] <- NA
      } else {
        z[i, j] <- zz
        start_array[, i, j] <- attr(zz, "free_pars")
        if (j == mle_idx2) {
          theta_start <- theta
        }
      }
    }
  }
  conf_region_list <- list(grid1 = grid1, grid2 = grid2, prof_loglik = z,
                           max_loglik = attr(object, "max_loglik"),
                           type = type, which_pars = which_pars,
                           name = attr(object, "name"))
  class(conf_region_list) <- "confreg"
  return(conf_region_list)
}

# ============================== conf_intervals ===============================

#' Confidence intervals
#'
#' Calculates confidence intervals for individual parameters
#'
#' @param object An object of class \code{"chandwich"} returned by
#'   \code{adjust_loglik}.
#' @param which_pars A vector specifying the (unfixed) parameters for which
#'   confidence intervals are required.  Can be either a numeric vector,
#'   specifying indices of the components of the \strong{full} parameter
#'   vector, or a character vector of parameter names, which must be a subset
#'   of those supplied in \code{par_names} in the call to
#'   \code{\link{adjust_loglik}} that produced \code{object}.
#'
#'   \code{which_pars} must not have any parameters in common with
#'   \code{attr(object, "fixed_pars")}.  \code{which_pars} must not contain
#'   all of the unfixed parameters, i.e. there is no point in profiling over
#'   all the unfixed parameters.
#' @param init A numeric vector of initial estimates of the values of the
#'   parameters that are not fixed and are not in \code{which_pars}.
#'   Should have length \code{attr(object, "p_current") - length(which_pars)}.
#'   If \code{init} is \code{NULL} or is of the wrong length then the
#'   relevant components from the MLE stored in \code{object} are used.
#' @param conf A numeric scalar in (0, 100). Confidence level for the
#'   intervals.
#' @param mult A numeric vector of length 1 or the same length as
#'   \code{which_pars}.
#'   The search for the profile loglikelihood-based confidence limits is
#'   conducted over the corresponding symmetric confidence intervals, extended
#'   by a factor of the corresponding component of \code{mult}.
#' @param num A numeric scalar.  The number of values at which to evaluate the
#'   profile loglikelihood either side of the MLE.
#' @param type A character scalar.  The argument \code{type} to the function
#'   returned by \code{\link{adjust_loglik}}, that is, the type of adjustment
#'   made to the independence loglikelihood function.
#' @param ... Further arguments to be passed to \code{\link[stats]{optim}}.
#'   These may include \code{gr}, \code{method}, \code{lower}, \code{upper}
#'   or \code{control}.
#' @return An object of class "confint", a list with components
#'     \item{conf}{The argument \code{conf}.}
#'     \item{cutoff}{A numeric scalar.  For values inside the
#'       confidence interval the profile loglikelihood lies above
#'       \code{cutoff}.}
#'     \item{parameter_vals, prof_loglik_vals}{\code{2 * num + 1} by
#'       \code{length{which_pars}} numeric matrices.
#'       Column i of \code{parameter_vals} contains the profiled values of
#'       parameter \code{which_par[i]}.  Column i of \code{prof_loglik_vals}
#'       contains the corresponding values of the profile loglikelihood.}
#'     \item{sym_CI, prof_CI}{\code{length(which_pars)}
#'       by 2 numeric matrices.  Row i of \code{sym_CI} (\code{prof_CI})
#'       contains the symmetric (profile loglikelihood-based) confidence
#'       intervals for parameter \code{which_pars[i]}.}
#'     \item{max_loglik}{The value of the adjusted loglikelihood
#'       at its maximum, stored in \code{object$max_loglik}.}
#'     \item{type}{The argument \code{type} supplied in the call
#'       to \code{conf_intervals}, i.e. the type of loglikelihood adjustment.}
#'     \item{which_pars}{The argument \code{which_pars}.}
#'     \item{name}{A character scalar. The name of the model,
#'       stored in \code{attr(object, "name")}.}
#'     \item{p_current}{The number of free parameters in the current model.}
#'     \item{fixed_pars, fixed_at}{\code{attr(object, "fixed_pars")} and
#'       \code{attr(object, "fixed_at")}, the arguments \code{fixed_pars} and
#'       \code{fixed_at} to \code{\link{adjust_loglik}}, if these were
#'       supplied.}
#' @seealso \code{\link{adjust_loglik}} to adjust a user-supplied
#'   loglikelihood function.
#' @seealso \code{\link{summary.chandwich}} for maximum likelihood estimates
#'   and unadjusted and adjusted standard errors.
#' @seealso \code{\link{plot.chandwich}} for plots of one-dimensional adjusted
#'   loglikelihoods.
#' @seealso \code{\link{conf_region}} for a confidence region for
#'   a pair of parameters.
#' @seealso \code{\link{compare_models}} to compare nested models using an
#'   (adjusted) likelihood ratio test.
#' @examples
#' # ------------------------- Binomial model, rats data ----------------------
#'
#' # Contributions to the independence loglikelihood
#' binom_loglik <- function(prob, data) {
#'   if (prob < 0 || prob > 1) {
#'     return(-Inf)
#'   }
#'   return(dbinom(data[, "y"], data[, "n"], prob, log = TRUE))
#' }
#' rat_res <- adjust_loglik(loglik = binom_loglik, data = rats, par_names = "p")
#'
#' # 95% likelihood-based confidence intervals, vertically adjusted
#' conf_intervals(rat_res)
#' \dontrun{
#' # Unadjusted
#' conf_intervals(rat_res, type = "none")
#' }
#'
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
#'                        par_names = par_names)
#'
#' # 95% likelihood-based confidence intervals, vertically adjusted
#' large_v <- conf_intervals(large, which_pars = c("xi[0]", "xi[1]"))
#' large_v
#' plot(large_v)
#' plot(large_v, which_par = "xi[1]")
#' \dontrun{
#' # Unadjusted
#' large_none <- conf_intervals(large, which_pars = c("xi[0]", "xi[1]"),
#'                              type = "none")
#' large_none
#' plot(large_v, large_none)
#' plot(large_v, large_none, which_par = "xi[1]")
#' }
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
#' conf_intervals(pois_quad)
#' @export
conf_intervals <- function(object, which_pars = NULL, init = NULL, conf = 95,
                     mult = 1.5, num = 10,
                     type = c("vertical", "cholesky", "spectral", "none"),
                     ...) {
  type <- match.arg(type)
  # Fixed parameters, values at which they are fixed and parameter names
  fixed_pars <- attr(object, "fixed_pars")
  fixed_at <- attr(object, "fixed_at")
  full_par_names <- attr(object, "full_par_names")
  # If which_pars is not supplied then set it to all (unfixed) parameters
  p <- attr(object, "p_full")
  if (is.null(which_pars)) {
    if (is.null(fixed_pars)) {
      which_pars <- 1:p
    } else {
      which_pars <- (1:p)[-fixed_pars]
    }
  }
  n_which_pars <- length(which_pars)
  n_mult <- length(mult)
  if (n_mult != 1 & n_mult != n_which_pars) {
    stop("mult must have length 1 or the same length as which_pars")
  } else if (n_mult == 1) {
    mult <- rep(mult, n_which_pars)
  }
  # If which_pars is a character vector
  if (is.character(which_pars)) {
    if (is.null(full_par_names)) {
      stop("which_pars can be character only if par_names is supplied")
    }
    if (!all(which_pars %in% full_par_names)) {
      stop("which_pars is not a subset of ", deparse(full_par_names))
    }
    temp <- which_pars
    which_pars <- which(full_par_names %in% which_pars)
    names(which_pars) <- temp
  } else {
    if (!is.null(full_par_names)) {
      names(which_pars) <- full_par_names[which_pars]
    }
  }
  if (any(which_pars %in% fixed_pars)) {
    stop("which_pars & attr(object,''fixed_pars'') have parameters in common")
  }
  # Symmetric confidence intervals
  if (is.null(fixed_pars)) {
    free_pars <- 1:p
  } else {
    free_pars <- (1:p)[-fixed_pars]
  }
  res_mle <- attr(object, "res_MLE")
  res_SE <- res_adjSE <- numeric(p)
  res_SE[free_pars] <- attr(object, "SE")
  res_adjSE[free_pars] <- attr(object, "adjSE")
  z_val <- stats::qnorm(1 - (1 - conf / 100) / 2)
  which_mle <- res_mle[which_pars]
  if (type == "none") {
    sym_lower <- which_mle - z_val * res_SE[which_pars]
    sym_upper <- which_mle + z_val * res_SE[which_pars]
    lower <- which_mle - mult * z_val * res_SE[which_pars]
    upper <- which_mle + mult * z_val * res_SE[which_pars]
  } else {
    sym_lower <- which_mle - z_val * res_adjSE[which_pars]
    sym_upper <- which_mle + z_val * res_adjSE[which_pars]
    lower <- which_mle - mult * z_val * res_adjSE[which_pars]
    upper <- which_mle + mult * z_val * res_adjSE[which_pars]
  }
  sym_CI <- cbind(sym_lower, sym_upper)
  colnames(sym_CI) <- c("lower", "upper")
  # Confidence intervals using profile loglikelihood
  prof_CI <- matrix(NA, nrow = n_which_pars, ncol = 2)
  colnames(prof_CI) <- c("lower", "upper")
  parameter_vals <- matrix(NA, nrow = 2 * num + 1, ncol = n_which_pars)
  prof_loglik_vals <- matrix(NA, nrow = 2 * num + 1, ncol = n_which_pars)
  if (!is.null(full_par_names)) {
    rownames(prof_CI) <- full_par_names[which_pars]
    colnames(parameter_vals) <- full_par_names[which_pars]
    colnames(prof_loglik_vals) <- full_par_names[which_pars]
  }
  # Insert the MLEs and the maximised loglikelihood values
  parameter_vals[num + 1, ] <- which_mle
  max_loglik <- attr(object, "max_loglik")
  prof_loglik_vals[num + 1, ] <- rep(max_loglik, n_which_pars)
  # 100% confidence interval: where the profile loglikelihood lies above cutoff
  cutoff <- max_loglik - stats::qchisq(conf  / 100, 1) / 2
  for (i in 1:length(which_pars)) {
    # MLE not including parameter being profiled and fixed parameters
    sol <- res_mle[-c(which_pars[i], fixed_pars)]
    # Lower tail ...
    par_low <- lower[i]
    par_vals <- seq(from = which_mle[i], to = par_low, length.out = num + 1)[-1]
    for (j in 1:num) {
      opt <- profile_loglik(object, prof_pars = which_pars[i],
                            prof_vals = par_vals[j],
                            init = sol, type = type, ...)
      sol <- attr(opt, "free_pars")
      parameter_vals[num - j + 1, i] <- par_vals[j]
      prof_loglik_vals[num - j + 1, i] <- opt
    }
    # Upper tail ...
    par_up <- upper[i]
    par_vals <- seq(from = which_mle[i], to = par_up, length.out = num + 1)[-1]
    for (j in 1:num) {
      opt <- profile_loglik(object, prof_pars = which_pars[i],
                            prof_vals = par_vals[j],
                            init = sol, type = type, ...)
      sol <- attr(opt, "free_pars")
      parameter_vals[num + j + 1, i] <- par_vals[j]
      prof_loglik_vals[num + j + 1, i] <- opt
    }
    # Use linear interpolation to estimate the confidence limits
    y <- prof_loglik_vals[, i]
    x <- parameter_vals[, i]
    # Find where the profile loglikelihood crosses cutoff
    temp <- diff(y - cutoff > 0)
    # Lower limit
    z <- which(temp == 1)
    low <- x[z] + (cutoff - y[z]) * (x[z + 1] - x[z]) / (y[z + 1] - y[z])
    # Upper limit
    z <- which(temp == -1)
    up <- x[z] + (cutoff - y[z]) * (x[z + 1] - x[z]) / (y[z + 1] - y[z])
    prof_CI[i, ] <- c(low, up)
  }
  conf_list <- list(conf = conf, cutoff = cutoff, parameter_vals = parameter_vals,
                    prof_loglik_vals = prof_loglik_vals, sym_CI = sym_CI,
                    prof_CI = prof_CI, max_loglik = attr(object, "max_loglik"),
                    type = type, which_pars = which_pars,
                    name = attr(object, "name"),
                    p_current = attr(object, "p_current"))
  if (!is.null(attr(object, "fixed_pars"))){
    conf_list <- c(conf_list, list(fixed_pars = attr(object, "fixed_pars"),
                                   fixed_at = attr(object, "fixed_at")))
  }
  class(conf_list) <- "confint"
  return(conf_list)
}

# ============================== profile_loglik ===============================

#' Profile loglikelihood
#'
#' Calculates the profile loglikelihood for a subset of the model parameters.
#' This function is provided primarily so that it can be called by
#' \code{\link{conf_intervals}} and \code{\link{conf_region}}.
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
#'   all of the unfixed parameters.
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
#'   loglikelihood.  The returned object has the attribute \code{"free_pars"}
#'   giving the optimal values of the parameters that remain after the
#'   parameters \code{prof_pars} and \code{attr(object, "fixed_pars")} have
#'   been removed from the full parameter vector.  If there are no such
#'   parameters, which happens if an attempt is made to profile over
#'   \emph{all} non-fixed parameters, then thiis attribute is not present and
#'   the value returned is calculated using the function \code{object}.
#' @seealso \code{\link{adjust_loglik}} to adjust a user-supplied
#'   loglikelihood function.
#' @seealso \code{\link{conf_intervals}} for confidence intervals for
#'   individual parameters.
#' @seealso \code{\link{conf_region}} for a confidence region for
#'   a pair of parameters.
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
#'                        par_names = par_names)
#'
#' # Profile loglikelihood for xi1, evaluated at xi1 = 0
#' profile_loglik(large, prof_pars = "xi[1]", prof_vals = 0)
#'
#' # Model with xi1 fixed at 0
#' medium <- adjust_loglik(larger = large, fixed_pars = "xi[1]")
#' # Profile loglikelihood for xi0, evaluated at xi0 = -0.1
#' profile_loglik(medium, prof_pars = "xi[0]", prof_vals = -0.1)
#' @export
profile_loglik <- function(object, prof_pars = NULL, prof_vals = NULL,
                           init = NULL, type = c("vertical", "cholesky",
                                                 "spectral", "none"), ...) {
  type <- match.arg(type)
  if (is.null(prof_pars)) {
    stop("prof_pars must be supplied")
  }
  fixed_pars <- attr(object, "fixed_pars")
  fixed_at <- attr(object, "fixed_at")
  full_par_names <- attr(object, "full_par_names")
  if (is.character(prof_pars)) {
    if (is.null(full_par_names)) {
      stop("prof_pars can be character only if par_names is supplied")
    }
    if (!all(prof_pars %in% full_par_names)) {
      stop("prof_pars is not a subset of ", deparse(full_par_names))
    }
    temp <- prof_pars
    prof_pars <- which(full_par_names %in% prof_pars)
    names(prof_pars) <- temp
  } else {
    if (!is.null(full_par_names)) {
      names(prof_pars) <- full_par_names[prof_pars]
    }
  }
  free_pars <- (1:attr(object, "p_full"))[-c(fixed_pars, prof_pars)]
  if (!is.null(full_par_names)) {
    names(free_pars) <- full_par_names[free_pars]
  }
  if (any(prof_pars %in% fixed_pars)) {
    stop("prof_pars & attr(object,''fixed_pars'') have parameters in common")
  }
  # The MLE (including any fixed parameters) of the full parameter vector
  full_mle <- attr(object, "res_MLE")
  if (is.null(prof_vals)) {
    prof_vals <- full_mle[prof_pars]
  }
  p_r <- length(free_pars)
  # If, after removing the parameters over which we wish to profile, there are
  # no free parameters, i.e. we are profiling over *all* unfixed parameters
  # then just use the adjusted loglikelihood given by type
  if (p_r == 0) {
    to_return <- do.call(object, list(prof_vals, type = type))
    return(to_return)
  }
  # Initial estimate: the MLE with fixed_pars set at the values in fixed_at
  if (is.null(init) || length(init) != p_r) {
    init <- full_mle[free_pars]
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
    pars[rank(c(prof_pars, free_pars))]<- c(prof_vals, x)
    return(-do.call(object, list(pars, type = type)))
  }
  # L-BFGS-B and Brent don't like Inf or NA or NaN
  if (optim_args$method == "L-BFGS-B" || optim_args$method == "Brent") {
    big_finite_val <- 10 ^ 10
    neg_prof_loglik <- function(x) {
      pars <- numeric(p)
      pars[rank(c(prof_pars, free_pars))]<- c(prof_vals, x)
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
  # Return the profile loglikelihood value (NOT negated)
  to_return <- -temp$value
  attr(to_return, "free_pars") <- temp$par
  return(to_return)
}
