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
