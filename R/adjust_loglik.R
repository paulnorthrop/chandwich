# ============================== adjust_loglik  ===============================
#
#' Loglikelihood adjustment
#'
#' Performs adjustments of a user-supplied independence loglikelihood for the
#' presence of cluster dependence, following
#' \href{http://dx.doi.org/10.1093/biomet/asm015}{Chandler and Bate (2007)}.
#' The user provides a function that returns observation-specifc
#' loglikelihood contributions and a vector that indicates cluster membership.
#'
#' @param loglik An R function.  Returns a vector of the
#'   loglikelihood contributions of individual observations.  The first
#'   argument must be the vector of model parameter(s). If any of the model
#'   parameters are out-of-bounds then \code{loglik} should return either
#'   \code{-Inf} or a vector with at least one element equal to \code{-Inf}.
#' @param ... Further arguments to be passed either to \code{loglik} or to
#'   \code{\link[stats]{optim}}.  The latter may include \code{gr},
#'   \code{method}, \code{lower}, \code{upper} or \code{control}.
#'   \code{hessian = TRUE} will be used regardless of any value supplied.
#'   The function \code{loglik} must \emph{not} have arguments with names
#'   that match any of these arguments to \code{\link[stats]{optim}}.
#' @param cluster A vector or factor indicating from which cluster the
#'   respective loglikelihood contributions from \code{loglik} originate.
#'   Must have the same length as the vector returned by \code{loglik}.
#'   By default each observation is its own cluster.
#' @param d A numeric scalar.  The number of model parameters.
#' @param init A numeric vector of initial values for use in the search for
#'   the MLE.  If \code{length(init)} is not equal to \code{d} then
#'   \code{length(init)} is taken as the number of model parameters.
#' @param par_names A character vector.  Names of the parameters.
#' @details Three adjustments to the independence loglikelihood are available.
#'   The `vertical' adjustment is described in Section 6 of
#'   Chandler and Bate (2007) and two `horizontal' adjustments are desribed in
#'   Sections 3.2 to 3.4 of Chandler and Bate (2007).
#'   See the descriptions of \code{type} and, for the
#'   horizontal adjustments, the descriptions of \code{C_cholesky} and
#'   \code{C_dilation}, in \strong{Value}.
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
#' @examples
#' binom_loglik <- function(prob, data) {
#'   if (prob < 0 || prob > 1) {
#'     return(-Inf)
#'   }
#'   return(dbinom(data[, "x"], data[, "size"], prob, log = TRUE))
#' }
#' prob <- 0.1
#' binom_data <- cbind(size = 10, x = rbinom(20, 10, prob))
#' cluster <- 1:nrow(binom_data)
#'
#' pjn <- adjust_loglik(loglik = binom_loglik, data = binom_data, cluster = cluster)
#' x <- seq(0.01, 0.99, by = 0.01)
#' y1 <- pjn(x, adjust = FALSE)
#' y2 <- pjn(x, type = "vertical")
#' y3 <- pjn(x, type = "cholesky")
#' y4 <- pjn(x, type = "dilation")
#' matplot(x, cbind(y1, y2, y3, y4), type = "l", lwd = 2)
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
#' @export
adjust_loglik <- function(loglik, ..., cluster = NULL, d = 1,
                          init = rep(0.1, d), par_names = NULL) {
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
  HI <- -temp$hessian
  for_grad <- list(func = neg_loglik, x = mle)
#  print("ZERO")
#  print(do.call(numDeriv::grad, for_grad))
  #
  # Find the estimated covariance matrix of the score vector ------------------
  #
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
