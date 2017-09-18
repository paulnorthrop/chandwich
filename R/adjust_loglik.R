# d, lower and upper, choose method (default based on d)
# transformation to [-Inf, Inf]?

# ============================== adjust_loglik  ===============================
#
#' Loglikelihood adjustment
#'
#' Chandler-Bate
#'
#' @param loglik An R function.  The first argument must be the vector of
#'   model parameter(s).  Returns a vector - use the length of this to check
#'   that cluster is the correct length!
#' @param ... Further arguments to be passed to logf.
#' @param cluster A numeric vector.  Must have length \code{length(data)} (if
#'   \code{data} is a vector) or \code{nrow(data)} (if \code{data} is a matrix.
#' @return A function that evaluates the adjusted loglikelihood.
#'    + ARGUMENTS horizontal1, horizontal2, vertical
#' @references Chandler, R. E. and Bate, S. (2007). Inference for clustered
#'   data using the independence loglikelihood. \emph{Biometrika},
#'   \strong{94}(1), 167-183. \url{http://dx.doi.org/10.1093/biomet/asm015}
#' @examples
#' binom_loglik <- function(prob, data) {
#'   return(dbinom(data[, "x"], data[, "size"], prob, log = TRUE))
#' }
#' prob <- 0.1
#' binom_data <- cbind(size = 10, x = rbinom(20, 10, prob))
#' binom_loglik(prob, binom_data)
#' cluster <- rep(1, nrow(binom_data))
#' adjust_loglik(loglik = binom_loglik, data = binom_data, cluster = cluster)
#' @export
adjust_loglik <- function(loglik, ..., cluster, d = 1, init = 0.1) {
  # Check that all the contributions to log-likelihood are finite at init
  check_vals <- loglik(init, ...)
  if (any(!is.finite(check_vals))) {
    stop("The loglikelihood is not finite at init")
  }
  # Check that cluster is a vector
  if (!is.vector(cluster)) {
    stop("cluster must be a vector")
  }
  # Check that cluster has the correct length
  n_cluster <- length(cluster)
  if (n_cluster != length(check_vals)) {
    stop("cluster must have the same length as the vector returned by loglik")
  }
  # Find the MLE
  neg_loglik <- function(x, ...) {
    return(-sum(loglik(x, ...)))
  }
  temp <- optim(init, neg_loglik, ..., lower = 0.01, upper = 0.99)
  return(temp)
}
