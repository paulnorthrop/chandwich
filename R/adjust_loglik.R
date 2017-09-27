# ============================== adjust_loglik  ===============================
#
#' Loglikelihood adjustment using the sandwich estimator
#'
#' Performs adjustments of a user-supplied independence loglikelihood for the
#' presence of cluster dependence, following
#' \href{http://dx.doi.org/10.1093/biomet/asm015}{Chandler and Bate (2007)}.
#' The user provides a function that returns observation-specifc
#' loglikelihood contributions and a vector that indicates cluster membership.
#'
#' @param loglik A function.  Returns a vector of the
#'   loglikelihood contributions of individual observations.  The first
#'   argument must be the vector of model parameter(s). If any of the model
#'   parameters are out-of-bounds then \code{loglik} should return either
#'   \code{-Inf} or a vector with at least one element equal to \code{-Inf}.
#' @param ... Further arguments to be passed either to \code{loglik}
#'   (and to \code{alg_deriv} and \code{alg_hess} if these are supplied) or
#'   to \code{\link[stats]{optim}}.  The latter may include \code{gr},
#'   \code{method}, \code{lower}, \code{upper} or \code{control}.
#'   \code{hessian = TRUE} will be used regardless of any value supplied.
#'   The function \code{loglik} must \emph{not} have arguments with names
#'   that match any of these arguments to \code{\link[stats]{optim}}.
#' @param cluster A vector or factor indicating from which cluster the
#'   respective loglikelihood contributions from \code{loglik} originate.
#'   Must have the same length as the vector returned by \code{loglik}.
#'   By default each observation is its own cluster.
#' @param d A numeric scalar.  The dimension of the parameter vector,
#'   i.e. the number of parameters in the model.
#' @param init A numeric vector of initial values for use in the search for
#'   the MLE.  If \code{length(init)} is not equal to \code{d} then
#'   \code{length(init)} is taken as the number of model parameters.
#' @param par_names A character vector.  Names of the parameters.
#' @param alg_deriv A function with the vector of model parameter(s) as its
#'   first argument.  Returns a \code{length(cluster)} by \code{d} numeric
#'   matrix. Column i contains the derivatives of each of the loglikelihood
#'   contributions in \code{loglik} with respect to model parameter i.
#' @param alg_hess A function with the vector of model parameter(s) as its
#'   first argument.  Returns a \code{d} by \code{d} numeric matrix equal to
#'   the Hessian of \code{loglik}, i.e. the matrix of second derivatives of
#'   the function \code{loglik}.
#' @details Three adjustments to the independence loglikelihood are available.
#'   The `vertical' adjustment is described in Section 6 of
#'   Chandler and Bate (2007) and two `horizontal' adjustments are described
#'   in Sections 3.2 to 3.4 of Chandler and Bate (2007).
#'   See the descriptions of \code{type} and, for the
#'   horizontal adjustments, the descriptions of \code{C_cholesky} and
#'   \code{C_dilation}, in \strong{Value}.
#'
#'   The adjustments involve first and second derviatives of the loglikelihood
#'   with respect to the model parameters.  These are estimated using
#'   \code{\link[numDeriv]{jacobian}} and \code{\link[stats]{optimHess}}
#'   unless \code{alg_deriv} and/or \code{alg_hess} are supplied.
#' @return A function of class \code{"chandwich"} to evaluate an adjusted
#'   loglikelihood, or the independence loglikelihood, at one or more sets
#'   of model parameters, with arguments
#'   \item{x}{A numeric vector or matrix of model parameter values.
#'     If \code{d = 1} this may be a numeric vector or a matrix with 1 column.
#'     If \code{d > 1} this may be a numeric vector of length \code{d}
#'     (one set of model parameters) or a numeric matrix with \code{d}
#'     columns (\code{ncol(x)} sets of model parameters).}
#'   \item{adjust}{A logical scalar.  Whether or not to adjust the
#'     independence loglikelihood.}
#'   \item{type}{A character scalar.  The type of adjustment to use if
#'     \code{adjust = TRUE}.  One of \code{"vertical"}, \code{"cholesky"} or
#'     \code{"dilation"}.}
#'   The function has (additional) attributes
#'   \item{d}{The number of models parameters.}
#'   \item{MLE}{The maximum likelihood estimate.}
#'   \item{SE}{The unadjusted standard errors.}
#'   \item{adjSE}{The adjusted standard errors.}
#'   \item{HI}{The Hessian of the independence loglikelihood.}
#'   \item{HA}{The Hessian of the adjusted loglikelihood.}
#'   \item{C_cholesky}{The matrix C in equation (14) of Chandler and Bate
#'     (2007), calculated using Cholesky decomposition.}
#'   \item{C_dilation}{The matrix C in equation (14) of Chandler and Bate
#'     (2007), calculated using spectral decomposition.}
#'   \item{par_names}{The argument \code{par_names}, if this was supplied.}
#' @references Chandler, R. E. and Bate, S. (2007). Inference for clustered
#'   data using the independence loglikelihood. \emph{Biometrika},
#'   \strong{94}(1), 167-183. \url{http://dx.doi.org/10.1093/biomet/asm015}
#' @seealso \code{\link{adjust_object}}: to adjust a fitted model object.
#' @seealso \code{\link{summary.chandwich}} for maximum likelihood estimates
#'   and unadjusted and adjusted standard errors.
#' @seealso \code{\link{plot.chandwich}} for one- and two- dimensional plots
#'   of of adjusted loglikelihoods.
#' @examples
#' # Binomial model, rats data ----------
#'
#' binom_loglik <- function(prob, data) {
#'   if (prob < 0 || prob > 1) {
#'     return(-Inf)
#'   }
#'   return(dbinom(data[, "y"], data[, "n"], prob, log = TRUE))
#' }
#' rat_res <- adjust_loglik(loglik = binom_loglik, data = rats)
#'
#' x <- seq(0.01, 0.99, by = 0.01)
#' y1 <- rat_res(x, adjust = FALSE)
#' y2 <- rat_res(x, type = "vertical")
#' y3 <- rat_res(x, type = "cholesky")
#' y4 <- rat_res(x, type = "dilation")
#' matplot(x, cbind(y1, y2, y3, y4), type = "l", lwd = 2)
#'
#' # Misspecified Poisson model for negative binomial data ----------
#' # ... following Section 5.1 of the
#' "Object-Oriented Computation of Sandwich Estimators" vignette of the
#' sandwich package
#' #https://cran.r-project.org/web/packages/sandwich/vignettes/sandwich-OOP.pdf
#'
#' set.seed(123)
#' x <- rnorm(250)
#' y <- rnbinom(250, mu = exp(1 + x), size = 1)
#' fm_pois <- glm(y ~ x + I(x^2), family = poisson)
#' adj_fn <- adjust_object(fm_pois)
#'
#' pois_glm_loglik <- function(pars, y, x) {
#'   log_mu <- pars[1] + pars[2] * x + pars[3] * x ^ 2
#'   return(dpois(y, lambda = exp(log_mu), log = TRUE))
#' }
#' pois_res <- adjust_loglik(pois_glm_loglik, y = y, x = x, d = 3)
#'
#' pois_alg_deriv <- function(pars, y, x) {
#'   mu <- exp(pars[1] + pars[2] * x + pars[3] * x ^ 2)
#'   return(cbind(y - mu, x * (y - mu), x ^2 * (y - mu)))
#' }
#'
#' pois_alg_hess <- function(pars, y, x) {
#'   mu <- exp(pars[1] + pars[2] * x + pars[3] * x ^ 2)
#'   alg_hess <- matrix(0, 3, 3)
#'   alg_hess[1, ] <- -c(sum(mu), sum(x * mu), sum(x ^ 2 * mu))
#'   alg_hess[2, ] <- -c(sum(x * mu), sum(x ^ 2 * mu), sum(x ^ 3 * mu))
#'   alg_hess[3, ] <- -c(sum(x ^ 2 * mu), sum(x ^ 3 * mu), sum(x ^ 4 * mu))
#'   return(alg_hess)
#' }
#'
#' pois_res <- adjust_loglik(pois_glm_loglik, y = y, x = x, d = 3,
#'                           alg_deriv = pois_alg_deriv, alg_hess = pois_alg_hess)
#'
#'
#' norm_loglik <- function(params, data) {
#'   mu <- params[1]
#'   sigma <- params[2]
#'   if (sigma <= 0) {
#'     return(-Inf)
#'   }
#'   return(dnorm(data, mean = mu, sd = sigma, log = TRUE))
#' }
#' mu <- 0
#' sigma <- 1
#' norm_data <- rnorm(2000, mean = mu, sd = sigma)
#' mu_sigma <- c(0, 1)
#' cluster <- 1:length(norm_data)
#' cluster <- rep(1:40, 50)
#'
#' pjn <- adjust_loglik(loglik = norm_loglik, data = norm_data, cluster = cluster,
#'               init = 0:1)
#'
#'
#' # GEV model, owtemps data ----------
#' # ... following Section 5.2 of Chandler and Bate (2007)
#'
#' if (requireNamespace("revdbayes", quietly = TRUE)) {
#'   gev_loglik <- function(pars, data) {
#'     o_pars <- pars[c(1, 3, 5)] + pars[c(2, 4, 6)]
#'     w_pars <- pars[c(1, 3, 5)] - pars[c(2, 4, 6)]
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
#' @export
adjust_loglik <- function(loglik, ..., cluster = NULL, d = 1,
                          init = rep(0.1, d), par_names = NULL,
                          alg_deriv = NULL, alg_hess = NULL) {
  #
  # Setup and checks -----------------------------------------------------------
  #
  # If d and length(init) don't agree then use length(init)
  if (d != length(init)) {
    d <- length(init)
  }
  # Extract from ... the arguments to be passed to stats::optim
  user_args <- list(...)
  optim_cond <- names(user_args) %in% methods::formalArgs(stats::optim)
  optim_args <- user_args[optim_cond]
  # The remaining arguments are to be passed loglik
  loglik_args <- user_args[!optim_cond]
  # Remove hessian, in case the user supplied it
  optim_args$hessian <- NULL
  # Check that all the contributions to loglikelihood are finite at init
  check_vals <- do.call(loglik, c(list(init), loglik_args))
  if (any(!is.finite(check_vals))) {
    stop("The loglikelihood is not finite at init")
  }
  # Number of terms in the loglikelihood
  n_loglik <- length(check_vals)
  # If cluster is not supplied then put observations in separate clusters
  # Otherwise, check that cluster is a vector of the correct length: n_loglik
  if (is.null(cluster)) {
    cluster <- 1:n_loglik
  } else {
    # Check that there are at least 2 clusters
    if (length(unique(cluster)) == 1) {
      stop("There must be more than one cluster")
    }
    # Check that cluster is a vector
    if (!is.vector(cluster) & !is.factor(cluster)) {
      stop("cluster must be a vector or a factor")
    }
    # Check that cluster has the correct length
    n_cluster <- length(cluster)
    if (n_cluster != length(check_vals)) {
      stop("cluster must have the same length as the vector returned by loglik")
    }
  }
  # Use "BFGS", unless the user has chosen the method or if d = 1 and they
  # have (inappropriately) chosen "Nelder-Mead" when d = 1
  if (is.null(optim_args$method)) {
    optim_args$method <- "BFGS"
  } else if (d == 1 & optim_args$method == "Nelder-Mead") {
    optim_args$method <- "BFGS"
  }
  # Define a function to minimise to find the MLE
  neg_loglik <- function(x) {
    loglik_vals <- do.call(loglik, c(list(x), loglik_args))
    return(-sum(loglik_vals))
  }
  # L-BFGS-B and Brent don't like Inf or NA or NaN
  if (optim_args$method == "L-BFGS-B" || optim_args$method == "Brent") {
    big_finite_val <- 10 ^ 10
    neg_loglik <- function(x) {
      loglik_vals <- do.call(loglik, c(list(x), loglik_args))
      check <- -sum(loglik_vals)
      if (!is.finite(check)) {
        check <- big_finite_val
      }
      return(check)
    }
  }
  for_optim <- c(list(par = init, fn = neg_loglik, hessian = TRUE), optim_args)
  #
  # Find the MLE and Hessian of the (negated) loglikelihood at the MLE -------
  #
  temp <- do.call(optim, for_optim)
  # Extract the MLE and the Hessian of independence loglikelihood at the MLE
  # Note the negation to change from Hessian of negated loglikelihood
  # to Hessian HI of loglikelihood
  mle <- temp$par
  max_loglik <- -temp$value
  if (is.null(alg_hess)) {
    HI <- -temp$hessian
  } else {
    HI <- alg_hess(mle, ...)
  }
#  for_grad <- list(func = neg_loglik, x = mle)
#  print("ZERO")
#  print(do.call(numDeriv::grad, for_grad))
  #
  # Find the estimated covariance matrix of the score vector ------------------
  #
  if (is.null(alg_deriv)) {
    # Function to aggregate the loglikelihood contributions within each cluster
    # [, 2] ensures that the correct *vector* is returned
    clus_loglik <- function(x, cluster) {
      loglik_vals <- do.call(loglik, c(list(x), loglik_args))
      return(aggregate(loglik_vals, list(cluster), sum)[, 2])
    }
    # Estimate the k x p matrix of derivatives of the k cluster-specific
    # loglikelihood contributions with respect to the p model parameters
    #
    for_jacobian <- list(func = clus_loglik, x = mle, cluster = cluster)
    U <- do.call(numDeriv::jacobian, for_jacobian)
  } else {
    U <- alg_deriv(mle, ...)
    U <- as.matrix(aggregate(U, list(cluster), sum)[, 2:(d+1)])
  }
  #
  # Unadjusted inverse Hessian and standard errors
  # [chol2inv(chol(X)) inverts X via its Cholesky decomposition]
  chol_minus_HI <- chol(-HI)
  HIinv <- -chol2inv(chol_minus_HI)
  SE <- sqrt(diag(-HIinv))
  # Adjusted Hessian and standard errors
  UHIinv <- U %*% HIinv
  HAinv <- -t(UHIinv) %*% UHIinv
  chol_minus_HAinv <- chol(-HAinv)
  HA <- -chol2inv(chol_minus_HAinv)
  adjSE <- sqrt(diag(-HAinv))
  # The following alternatives give the same answer ...
  # Estimate covariance of score using equation (7) of Chandler and Bate (2007)
  #   V <- t(U) %*% U
  #   Vinv <- chol2inv(chol(V))
  #   HA <- -solve(HIinv %*% V %*% HIinv)
  #   HA <- -HI %*% Vinv %*% HI
  #
  # Calculate the matrix C used in the horizontal adjustment described in
  # Section 3.2 of CB2007.  We do this in two ways, using Cholesky and
  # spectral decompositions, respectively (see Section 3.4 of CB2007).
  #
  # Note that we decompose the *negated* Hessians of the loglikelihood
  # because the Hessians are not positive definite.
  MI <- chol_minus_HI
  MA <- chol(-HA)
  C_cholesky <- solve(MI) %*% MA
  z <- eigen(-HI, symmetric = TRUE)
  # We need nrow and ncol arguments to diag() so that the d = 1 case is correct
  MI <- z$vectors %*% diag(sqrt(z$values), d, d) %*% t(z$vectors)
  z <- eigen(-HA, symmetric = TRUE)
  MA <- z$vectors %*% diag(sqrt(z$values), d, d) %*% t(z$vectors)
  C_dilation <- solve(MI) %*% MA
  # Return a function to calculate the adjusted loglikelihood
  # If x = mle then we return the maximum of the independence loglikelihood
  adjust_loglik_fn <- function(x, adjust = TRUE,
                               type = c("vertical", "cholesky", "dilation")) {
    type <- match.arg(type)
    x <- as.matrix(x)
    if (d > 1 & ncol(x) == 1) {
      x <- t(x)
    }
    if (ncol(x) != d) {
      stop("x does not have the correct dimensions")
    }
    if (!adjust) {
      fn <- function(x) {
        if (identical(x, mle)) {
          return(max_loglik)
        }
        loglik_vals <- do.call(loglik, c(list(x), loglik_args))
        return(sum(loglik_vals))
      }
      return(apply(x, 1, fn))
    }
    if (type == "vertical") {
      fn <- function(x) {
        if (identical(x, mle)) {
          return(max_loglik)
        }
        loglik_vals <- do.call(loglik, c(list(x), loglik_args))
        ind_loglik <- sum(loglik_vals)
        snum <- t(x - mle) %*% HA %*% (x - mle)
        sden <- t(x - mle) %*% HI %*% (x - mle)
        s <- snum / sden
        return(max_loglik + s * (ind_loglik - max_loglik))
      }
      return(apply(x, 1, fn))
    } else {
      fn <- function(x) {
        if (identical(x, mle)) {
          return(max_loglik)
        }
        C <- ifelse(type == "cholesky", C_cholesky, C_dilation)
        x_star <- mle + (C %*% (x - mle))
        loglik_vals <- do.call(loglik, c(list(x_star), loglik_args))
        return(sum(loglik_vals))
      }
      return(apply(x, 1, fn))
    }
  }
  attr(adjust_loglik_fn, "d") <- d
  attr(adjust_loglik_fn, "MLE") <- mle
  attr(adjust_loglik_fn, "SE") <- SE
  attr(adjust_loglik_fn, "adjSE") <- adjSE
  attr(adjust_loglik_fn, "HI") <- HI
  attr(adjust_loglik_fn, "HA") <- HA
  attr(adjust_loglik_fn, "C_cholesky") <- C_cholesky
  attr(adjust_loglik_fn, "C_dilation") <- C_dilation
  attr(adjust_loglik_fn, "par_names") <- par_names
  class(adjust_loglik_fn) <- "chandwich"
  return(adjust_loglik_fn)
}
