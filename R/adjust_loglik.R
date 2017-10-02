# ============================== adjust_loglik  ===============================
#
#' Loglikelihood adjustment using the sandwich estimator
#'
#' Performs adjustments of a user-supplied independence loglikelihood for the
#' presence of cluster dependence, following
#' \href{http://dx.doi.org/10.1093/biomet/asm015}{Chandler and Bate (2007)}.
#' The user provides a function that returns observation-specifc
#' loglikelihood contributions and a vector that indicates cluster membership.
#' The loglikelihood of a sub-model can be adjusted by fixing a set of
#' parameters at particular values.
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
#'   In the call to \code{\link[stats]{optim}}, \code{hessian = TRUE}
#'   will be used regardless of any value supplied.
#'   The function \code{loglik} must \emph{not} have arguments with names
#'   that match any of these arguments to \code{\link[stats]{optim}}.
#' @param cluster A vector or factor indicating from which cluster the
#'   respective loglikelihood contributions from \code{loglik} originate.
#'   Must have the same length as the vector returned by \code{loglik}.
#'   By default each observation is its own cluster.
#' @param p A numeric scalar.  The dimension of the \emph{full} parameter
#'   vector, i.e. the number of parameters in the full model.
#' @param init A numeric vector of initial values.  Must have length equal
#'   to the number of parameters in the full model.  If \code{init} is
#'   supplied then \code{p} is set to \code{length(init)}.
#'   If \code{fixed_pars} is not \code{NULL} then \code{init[-fixed_pars]}
#'   is used in the search for the MLE.
#'   If \code{init} is not supplied then \code{rep(0.1, p)} is used.
#' @param par_names A character vector.  Names of the \code{p} parameters
#'   in the full model.
#' @param fixed_pars A numeric vector.  Indices of the components of the
#'   parameter vector that are restricted to be equal to the value(s) in
#'   \code{fixed_at}.
#' @param fixed_at A numeric vector of length 1 or \code{length(fixed_pars)}.
#'   If \code{length(fixed_at) = 1} then the components \code{fixed_pars}
#'   of the parameter vector are all fixed at \code{fixed_at}.
#'   If \code{length(fixed_at) = length(fixed_pars)} then the component
#'   \code{fixed_pars[i]} is fixed at \code{fixed_at[i]} for each \code{i}.
#' @param larger (Only relevant if \code{fixed_pars} is not \code{NULL}.)
#'   An object of class \code{"chandwich"} returned by \code{adjust_loglik},
#'   corresponding to a model in which the smaller model implied by
#'   \code{fixed_pars} is nested.  If \code{larger} is supplied then
#'   all the arguments to \code{adjust_loglik} apart from
#'   \code{fixed_pars} and \code{fixed_at} are extracted from \code{larger}.
#'   If \code{init} is not supplied in the current call to
#'   \code{adjust_loglik} then \code{init} is set to
#'   \code{attr(larger, "MLE")}, with the elements in \code{fixed_pars}
#'   set to \code{fixed_at}.
#' @param alg_deriv A function with the vector of model parameter(s) as its
#'   first argument.  Returns a \code{length(cluster)} by \code{p} numeric
#'   matrix. Column i contains the derivatives of each of the loglikelihood
#'   contributions in \code{loglik} with respect to model parameter i.
#' @param alg_hess A function with the vector of model parameter(s) as its
#'   first argument.  Returns a \code{p} by \code{p} numeric matrix equal to
#'   the Hessian of \code{loglik}, i.e. the matrix of second derivatives of
#'   the function \code{loglik}.
#' @details Three adjustments to the independence loglikelihood described in
#'   Chandler and Bate (2007) are available.  The `vertical' adjustment is
#'   described in Section 6 and two `horizontal' adjustments are described
#'   in Sections 3.2 to 3.4.  See the descriptions of \code{type} and, for the
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
#'   \item{x}{A numeric vector or matrix giving values of the \code{p_current}
#'     (see below) parameters in the model to which the returned adjusted
#'     loglikelihood applies.
#'     If \code{p = 1} this may be a numeric vector or a matrix with 1 column.
#'     If \code{p > 1} this may be a numeric vector of length \code{p}
#'     (one set of model parameters) or a numeric matrix with \code{p}
#'     columns (\code{ncol(x)} sets of model parameters).}
#'   \item{type}{A character scalar.  The type of adjustment to use.
#'     One of \code{"vertical"}, \code{"cholesky"}, \code{"dilation"} or
#'     \code{"none"}.}  The latter results in the evaluation of the
#'     (unadjusted) independence loglikelihood.
#'   The function has (additional) attributes
#'   \item{p_full, p_current}{The number of parameters in the full model and
#'     current models, respectively.}
#'   \item{MLE}{The maximum likelihood estimate.}
#'   \item{res_MLE}{The maximum likelihood estimate, including any parameters
#'     with fixed values. Equal to MLE if \code{fixed_pars} was \code{NULL}.}
#'   \item{SE, adjSE}{The unadjusted and adjusted standard errors, respectively.}
#'   \item{HI, HA}{The Hessians of the independence and adjusted loglikelihood,
#'     respectively.}
#'   \item{C_cholesky, C_dilation}{The matrix C in equation (14) of Chandler and
#'     Bate (2007), calculated using Cholesky decomposition and spectral
#'     decomposition, respectively.}
#'   \item{full_par_names, par_names}{The names of the parameters in the full
#'     and current models, respectively, if these were supplied in
#'     this call or a previous call.}
#'   \item{max_loglik}{The common maximised value of the independence and
#'     adjusted loglikelihoods.}
#'   \item{loglik, cluster}{The arguments \code{loglik} and \code{cluster}
#'     supplied in this call, or a previous call.}
#'   \item{loglik_args}{A list containing the further arguments passed to
#'     \code{loglik} via ... in this call, or a previous call.}
#'   If \code{fixed_pars} is not \code{NULL} then there are further attributes
#'   \item{fixed_pars}{The argument \code{fixed_pars}.}
#'   \item{fixed_at}{The argument \code{fixed_at}.}
#'   If \code{alg_deriv} and/or \code{alg_hess} were supplied then these are
#'   returned as further attributes.
#' @references Chandler, R. E. and Bate, S. (2007). Inference for clustered
#'   data using the independence loglikelihood. \emph{Biometrika},
#'   \strong{94}(1), 167-183. \url{http://dx.doi.org/10.1093/biomet/asm015}
#' @seealso \code{\link{summary.chandwich}} for maximum likelihood estimates
#'   and unadjusted and adjusted standard errors.
#' @seealso \code{\link{plot.chandwich}} for one- and two- dimensional plots
#'   of of adjusted loglikelihoods.
#' @seealso \code{\link{compare_models}} to compare nested models using an
#'   (adjusted) likelihood ratio test.
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
#' y1 <- rat_res(x, type = "none")
#' y2 <- rat_res(x, type = "vertical")
#' y3 <- rat_res(x, type = "cholesky")
#' y4 <- rat_res(x, type = "dilation")
#' matplot(x, cbind(y1, y2, y3, y4), type = "l", lwd = 2)
#'
#' # Misspecified Poisson model for negative binomial data ----------
#' # ... following Section 5.1 of the
#' # "Object-Oriented Computation of Sandwich Estimators" vignette of the
#' # sandwich package
#' # https://cran.r-project.org/web/packages/sandwich/vignettes/sandwich-OOP.pdf
#'
#' set.seed(123)
#' x <- rnorm(250)
#' y <- rnbinom(250, mu = exp(1 + x), size = 1)
#' fm_pois <- stats::glm(y ~ x + I(x^2), family = poisson)
#'
#' pois_glm_loglik <- function(pars, y, x) {
#'   log_mu <- pars[1] + pars[2] * x + pars[3] * x ^ 2
#'   return(dpois(y, lambda = exp(log_mu), log = TRUE))
#' }
#' pois_res <- adjust_loglik(pois_glm_loglik, y = y, x = x, p = 3)
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
#' pois_res <- adjust_loglik(pois_glm_loglik, y = y, x = x, p = 3,
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
#'   gev_loglik <- function(pars, data) {
#'     o_pars <- pars[c(1, 3, 5)] + pars[c(2, 4, 6)]
#'     w_pars <- pars[c(1, 3, 5)] - pars[c(2, 4, 6)]
#'     if (o_pars[2] <= 0 | w_pars[2] <= 0) return(-Inf)
#'     o_loglik <- gev_dens(data[, "Oxford"], o_pars[1], o_pars[2], o_pars[3],
#'                          log = TRUE)
#'     w_loglik <- gev_dens(data[, "Worthing"], w_pars[1], w_pars[2], w_pars[3],
#'                          log = TRUE)
#'     return(o_loglik + w_loglik)
#'   }
#'
#' # Initial estimates (method of moments for the Gumbel case)
#' sigma <- as.numeric(sqrt(6 * diag(stats::var(owtemps))) / pi)
#' mu <- as.numeric(colMeans(owtemps) - 0.57722 * sigma)
#' init <- c(mean(mu), -diff(mu) / 2, mean(sigma), -diff(sigma) / 2, 0, 0)
#'
#' # Perform the log-likelihood adjustment
#' large <- adjust_loglik(gev_loglik, data = owtemps, init = init,
#'         par_names = c("mu0", "mu1", "sigma0", "sigma1", "xi0", "xi1"))
#' # Rows 1, 3 and 4 of Table 2 of Chandler and Bate (2007)
#' round(attr(large, "MLE"), 4)
#' round(attr(large, "SE"), 4)
#' round(attr(large, "adjSE"), 4)
#'
#' # Easier way to fix parameters
#' medium <- adjust_loglik(larger = large, fixed_pars = 6)
#' small1 <- adjust_loglik(larger = large, fixed_pars = c(4, 6))
#' # or
#' small2 <- adjust_loglik(larger = medium, fixed_pars = c(4, 6))
#'
#' # (Slightly) harder way
#' medium <- adjust_loglik(gev_loglik, data = owtemps, init = init, fixed_pars = 6,
#'           par_names = c("mu0", "mu1", "sigma0", "sigma1", "xi0", "xi1"))
#' small <- adjust_loglik(gev_loglik, data = owtemps, init = init, fixed_pars = c(4,6),
#'          par_names = c("mu0", "mu1", "sigma0", "sigma1", "xi0", "xi1"))
#' @export
adjust_loglik <- function(loglik = NULL, ..., cluster = NULL, p = 1,
                          init = NULL, par_names = NULL,
                          fixed_pars = NULL, fixed_at = 0, larger = NULL,
                          alg_deriv = NULL, alg_hess = NULL) {
  #
  # Setup and checks -----------------------------------------------------------
  #
  if (is.null(loglik) & is.null(larger)) {
    stop("If loglik is NULL then larger (and fixed_pars) must be supplied")
  }
  # If larger is NULL then extract all information from the supplied arguments.
  # Otherwise, use the information contained in larger and only allow the user
  # to override the initial estimates init.
  if (is.null(larger)) {
    got_loglik_args <- FALSE
    if (!is.null(init)) {
      p <- length(init)
    } else {
      init <- rep(0.1, p)
    }
    if (!is.null(fixed_pars)) {
      init <- init[-fixed_pars]
    }
    # If par_names is supplied then force it to have length p
    if (!is.null(par_names)) {
      full_par_names <- par_names <- rep_len(par_names, p)
      if (!is.null(fixed_pars)) {
        par_names <- par_names[-fixed_pars]
      }
    } else {
      full_par_names <- NULL
    }
  } else {
    if (is.null(fixed_pars)) {
      warning("larger is only relevant if fixed_pars is supplied")
      got_loglik_args <- FALSE
      full_par_names <- NULL
    } else {
      if (!inherits(larger, "chandwich")) {
        stop("larger must be a \"chandwich\" object")
      }
      loglik <- attr(larger, "loglik")
      loglik_args <- attr(larger, "loglik_args")
      got_loglik_args <- TRUE
      cluster <- attr(larger, "cluster")
      p <- attr(larger, "p_full")
      if (is.null(init)) {
        init <- attr(larger, "res_MLE")[-fixed_pars]
      } else {
        init <- init[-fixed_pars]
      }
      full_par_names <- attr(larger, "full_par_names")
      par_names <- attr(larger, "full_par_names")[-fixed_pars]
      alg_deriv <- attr(larger, "alg_deriv")
      alg_hess <- attr(larger, "alg_hess")
    }
  }
  # Extract from ... the arguments to be passed to stats::optim
  user_args <- list(...)
  optim_cond <- names(user_args) %in% methods::formalArgs(stats::optim)
  optim_args <- user_args[optim_cond]
  # The remaining arguments are to be passed loglik
  # Only extract these if they were not extracted from larger earlier
  if (!got_loglik_args) {
    loglik_args <- user_args[!optim_cond]
  }
  # Remove hessian, in case the user supplied it
  optim_args$hessian <- NULL
  # Set the number of parameters and the initial estimates
  if (is.null(fixed_pars)) {
    n_pars <- p
    # Check that all the contributions to loglikelihood are finite at init
    check_vals <- do.call(loglik, c(list(init), loglik_args))
  } else {
    qq <- length(fixed_pars)
    n_pars <- p - qq
    if (qq >= p) {
      stop("length(fixed_pars) must be smaller than p")
    }
    if (!(length(fixed_at) %in% c(1, qq))) {
      stop("the lengths of 'fixed_pars' and 'fixed_at' are not compatible")
    }
    fixed_at <- rep_len(fixed_at, qq)
    free_pars <- (1:p)[-fixed_pars]
#    par_names <- par_names[free_pars]
    # Check that all the contributions to the loglikelihood are finite at init
    pars <- numeric(p)
    pars[fixed_pars] <- fixed_at
    pars[free_pars] <- init
    check_vals <- do.call(loglik, c(list(pars), loglik_args))
  }
  if (any(!is.finite(check_vals))) {
    stop("The loglikelihood is not finite at init")
  }
  # Number of terms in the loglikelihood
  n_loglik <- length(check_vals)
  if (n_loglik == 1) {
    stop("There must be more than one cluster")
  }
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
  # Use "BFGS", unless the user has chosen the method or if p = 1 and they
  # have (inappropriately) chosen "Nelder-Mead" when p = 1
  if (is.null(optim_args$method)) {
    optim_args$method <- "BFGS"
  } else if (p == 1 & optim_args$method == "Nelder-Mead") {
    optim_args$method <- "BFGS"
  }
  if (is.null(fixed_pars)) {
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
  } else {
    neg_loglik <- function(x) {
      pars <- numeric(p)
      pars[fixed_pars] <- fixed_at
      pars[free_pars] <- x
      loglik_vals <- do.call(loglik, c(list(pars), loglik_args))
      return(-sum(loglik_vals))
    }
    # L-BFGS-B and Brent don't like Inf or NA or NaN
    if (optim_args$method == "L-BFGS-B" || optim_args$method == "Brent") {
      big_finite_val <- 10 ^ 10
      neg_loglik <- function(x) {
        pars <- numeric(p)
        pars[fixed_pars] <- fixed_at
        pars[free_pars] <- x
        loglik_vals <- do.call(loglik, c(list(pars), loglik_args))
        check <- -sum(loglik_vals)
        if (!is.finite(check)) {
          check <- big_finite_val
        }
        return(check)
      }
    }
  }
  for_optim <- c(list(par = init, fn = neg_loglik, hessian = TRUE), optim_args)
  #
  # Find the MLE and Hessian of the (negated) loglikelihood at the MLE -------
  #
  temp <- do.call(stats::optim, for_optim)
  # Extract the MLE and the Hessian of independence loglikelihood at the MLE
  # Note the negation to change from Hessian of negated loglikelihood
  # to Hessian HI of loglikelihood
  mle <- temp$par
  if (!is.null(fixed_pars)) {
    res_mle <- numeric(p)
    res_mle[fixed_pars] <- fixed_at
    res_mle[free_pars] <- mle
  }
  max_loglik <- -temp$value
  if (is.null(alg_hess)) {
    HI <- -temp$hessian
  } else {
    if (is.null(fixed_pars)) {
      HI <- do.call(alg_hess, c(list(mle), loglik_args))
    } else {
      HI <- do.call(alg_hess, c(list(res_mle), loglik_args))
      HI <- HI[-fixed_pars, -fixed_pars]
    }
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
    if (is.null(fixed_pars)) {
        clus_loglik <- function(x, cluster) {
          loglik_vals <- do.call(loglik, c(list(x), loglik_args))
          return(stats::aggregate(loglik_vals, list(cluster), sum)[, 2])
        }
    } else {
      clus_loglik <- function(x, cluster) {
        pars <- numeric(p)
        pars[fixed_pars] <- fixed_at
        pars[free_pars] <- x
        loglik_vals <- do.call(loglik, c(list(pars), loglik_args))
        return(stats::aggregate(loglik_vals, list(cluster), sum)[, 2])
      }
    }
    # Estimate the k x p matrix of derivatives of the k cluster-specific
    # loglikelihood contributions with respect to the p model parameters
    #
    for_jacobian <- list(func = clus_loglik, x = mle, cluster = cluster)
    U <- do.call(numDeriv::jacobian, for_jacobian)
  } else {
    if (is.null(fixed_pars)) {
      U <- do.call(alg_deriv, c(list(mle), loglik_args))
      U <- as.matrix(stats::aggregate(U, list(cluster), sum)[, 2:(p + 1)])
    } else {
      U <- do.call(alg_deriv, c(list(res_mle), loglik_args))
      U <- as.matrix(stats::aggregate(U, list(cluster), sum)[, 2:(p + 1)])
      U <- U[, -fixed_pars]
    }
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
  MI <- z$vectors %*% diag(sqrt(z$values), n_pars, n_pars) %*% t(z$vectors)
  z <- eigen(-HA, symmetric = TRUE)
  MA <- z$vectors %*% diag(sqrt(z$values), n_pars, n_pars) %*% t(z$vectors)
  C_dilation <- solve(MI) %*% MA
  # If some parameters are fixed then modify the input loglik so that it
  # accepts an argument of length length(free_pars)
  if (!is.null(fixed_pars)) {
    ret_loglik <- function(x) {
      pars <- numeric(p)
      pars[fixed_pars] <- fixed_at
      pars[free_pars] <- x
      return(do.call(loglik, c(list(pars), loglik_args)))
    }
  } else {
    ret_loglik <- function(x) {
      return(do.call(loglik, c(list(x), loglik_args)))
    }
  }
  # Return a function to calculate the adjusted loglikelihood
  # If x = mle then we return the maximum of the independence loglikelihood
  adjust_loglik_fn <- function(x, type = c("vertical", "cholesky", "dilation",
                                           "none")) {
    type <- match.arg(type)
    x <- as.matrix(x)
    if (n_pars > 1 & ncol(x) == 1) {
      x <- t(x)
    }
    if (ncol(x) != n_pars) {
      stop("x does not have the correct dimensions")
    }
    if (type == "vertical") {
      fn <- function(x) {
        if (identical(x, mle)) {
          return(max_loglik)
        }
        loglik_vals <- do.call(ret_loglik, list(x))
        ind_loglik <- sum(loglik_vals)
        snum <- t(x - mle) %*% HA %*% (x - mle)
        sden <- t(x - mle) %*% HI %*% (x - mle)
        s <- snum / sden
        return(max_loglik + s * (ind_loglik - max_loglik))
      }
      return(apply(x, 1, fn))
    } else if (type %in% c("cholesky", "dilation")) {
      fn <- function(x) {
        if (identical(x, mle)) {
          return(max_loglik)
        }
        C <- ifelse(type == "cholesky", C_cholesky, C_dilation)
        x_star <- mle + (C %*% (x - mle))
        loglik_vals <- do.call(ret_loglik, list(x_star))
        return(sum(loglik_vals))
      }
      return(apply(x, 1, fn))
    } else {
      fn <- function(x) {
        if (identical(x, mle)) {
          return(max_loglik)
        }
        loglik_vals <- do.call(ret_loglik, list(x))
        return(sum(loglik_vals))
      }
      return(apply(x, 1, fn))
    }
  }
  attr(adjust_loglik_fn, "p_full") <- p
  attr(adjust_loglik_fn, "p_current") <- n_pars
  if (!is.null(fixed_pars)) {
    attr(adjust_loglik_fn, "fixed_pars") <- fixed_pars
    attr(adjust_loglik_fn, "fixed_at") <- fixed_at
    attr(adjust_loglik_fn, "res_MLE") <- res_mle
  } else {
    attr(adjust_loglik_fn, "res_MLE") <- mle
  }
  attr(adjust_loglik_fn, "alg_deriv") <- alg_deriv
  attr(adjust_loglik_fn, "alg_hess") <- alg_hess
  attr(adjust_loglik_fn, "MLE") <- mle
  attr(adjust_loglik_fn, "SE") <- SE
  attr(adjust_loglik_fn, "adjSE") <- adjSE
  attr(adjust_loglik_fn, "HI") <- HI
  attr(adjust_loglik_fn, "HA") <- HA
  attr(adjust_loglik_fn, "C_cholesky") <- C_cholesky
  attr(adjust_loglik_fn, "C_dilation") <- C_dilation
  attr(adjust_loglik_fn, "par_names") <- par_names
  attr(adjust_loglik_fn, "full_par_names") <- full_par_names
  attr(adjust_loglik_fn, "loglik") <- loglik
  attr(adjust_loglik_fn, "cluster") <- cluster
  attr(adjust_loglik_fn, "max_loglik") <- max_loglik
  attr(adjust_loglik_fn, "loglik_args") <- loglik_args
  class(adjust_loglik_fn) <- "chandwich"
  return(adjust_loglik_fn)
}
