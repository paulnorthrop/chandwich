# ============================== plot.chandwich ===============================

#' Plot diagnostics a chandwich object
#'
#' \code{plot} method for class "chandwich".
#'
#' @param x an object of class "chandwich", a result of a call to
#'   \code{\link{adjust_loglik}}.
#' @param y Not used.
#' @param lower,upper Numeric vectors specifying the lower and upper limits
#'   on the parameter values to appear in the plot.  If either \code{lower}
#'   or \code{upper} are not provided then the MLE minus (for \code{lower})
#'   or plus (for \code{upper}) three (adjusted) standard errors is used.
#' @param type An integer vector, a subset of the numbers \code{1:4}.
#'   Indicates which loglikelihoods to plot: \code{1} for \code{"vertical"}
#'   adjustment; \code{2} for \code{"cholesky"} (horizontal adjustment);
#'   \code{3} for \code{"dilation"} (horizontal adjustment); \code{4}
#'   for no adjustment, i.e. based on the independence loglikelihood.
#' @param ... Additional arguments passed on to ...
#' @return A list containing the graphical parameters using in producing the
#'   plot including any arguments supplied via ... is returned (invisibly).
#' @examples
#' binom_loglik <- function(prob, data) {
#'   if (prob < 0 || prob > 1) {
#'     return(-Inf)
#'   }
#'   return(dbinom(data[, "y"], data[, "n"], prob, log = TRUE))
#' }
#' cluster <- 1:nrow(rats)
#' rat_res <- adjust_loglik(loglik = binom_loglik, data = rats, cluster = cluster)
#' plot(rat_res)
#' plot(rat_res, type = 1:4)
#' plot(rat_res, type = 1:4, lower = 0, upper = 1)
#' @seealso \code{\link{adjust_loglik}}.
#' @seealso \code{\link{summary.chandwich}} for maximum likelihood estimates
#'   and unadjusted and adjusted standard errors.
#' @export
plot.chandwich <- function(x, y, ..., lower = NULL, upper = NULL,
                           type = 1) {
  if (!inherits(x, "chandwich")) {
    stop("use only with \"chandwich\" objects")
  }
  n_pars <- attr(x, "p_current")
  # Single parameter model
  if (n_pars == 1) {
    if (is.null(lower)) {
      lower <- attr(x, "MLE") - 3 * attr(x, "adjSE")
      upper <- attr(x, "MLE") + 3 * attr(x, "adjSE")
    }
    x_vals <- seq(lower, upper, len = 100)
    y <- NULL
    if (any(type == 1)) {
      y <- cbind(y, x(x_vals, type = "vertical"))
    }
    if (any(type == 2)) {
      y <- cbind(y, x(x_vals, type = "cholesky"))
    }
    if (any(type == 3)) {
      y <- cbind(y, x(x_vals, type = "dilation"))
    }
    if (any(type == 4)) {
      y <- cbind(y, x(x_vals, adjust = FALSE))
    }
    graphics::matplot(x_vals, y, type = "l", ...)
  }
  return(invisible())
}

# ============================ summary.chandwich =============================

#' Summarizing adjusted loglikelihoods
#'
#' \code{summary} method for class "chandwich"
#'
#' @param object an object of class "chandwich", a result of a call to
#'   \code{adjust_loglik}.
#' @param digits A integer. Used for number formatting with
#'   \code{\link{signif}}.
#' @param ... Additional optional arguments. At present no optional
#'   arguments are used.
#' @return Returns a numeric matrix with 3 or 4 columns and \code{object$d}
#'   rows, where \code{object$d} is the number of parameters in the model
#'   defined by the argument \code{loglik} in the call to
#'   \code{\link{adjust_loglik}}.
#'   The columns contain: the parameter names (if these were supplied using
#'   \code{par_names} in the call to \code{\link{adjust_loglik}}),
#'   the MLE, unadjusted standard errors (SE) and adjusted standard errors
#'   (adjSE).
#' @examples
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
#' summary(pjn)
#' @seealso \code{\link{adjust_loglik}}.
#' @export
summary.chandwich <- function(object, digits = max(3, getOption("digits")-3),
                              ...) {
  if (!inherits(object, "chandwich")) {
    stop("use only with \"chandwich\" objects")
  }
  mle <- signif(attr(object, "MLE"), digits = digits)
  SE <- signif(attr(object, "SE"), digits = digits)
  adjSE <- signif(attr(object, "adjSE"), digits = digits)
  res <- cbind(mle, SE, adjSE)
  column_names <- c("MLE", "SE", "adjSE")
  if (!is.null(attr(object, "par_names"))) {
    res <- cbind(attr(object, "par_names"), res)
    column_names <- c("parameter", column_names)
  }
  colnames(res) <- column_names
#  rownames(res) <- 1:length(object$v_vec)
  return(res)
}
