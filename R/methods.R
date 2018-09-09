# ============================== coef.chandwich ==============================

#' Extract model coefficients method for objects of class "chandwich"
#'
#' \code{coef} method for class "chandwich".
#'
#' @param object an object of class "chandwich", a result of a call to
#'   \code{\link{adjust_loglik}}.
#' @param complete A logical scalar.  If \code{complete = TRUE} then the
#'   full vector of parameter estimates is returned, including any fixed
#'   using \code{fixed_pars} and \code{fixed_at} in the call to
#'   \code{\link{adjust_loglik}}.  Otherwise, only the vector of estimates
#'   of free parameters is returned.
#' @param ... Additional optional arguments. At present no optional
#'   arguments are used.
#' @details The full vector of estimate is taken from
#' \code{attributes(object)$res_MLE} and the reduced vector from
#' \code{attributes(object)$MLE}.
#' @return The argument \code{x}, invisibly, as for all \code{\link{print}}
#'   methods.
#' @seealso \code{\link{summary.chandwich}}: \code{summary} method for
#'   class "chandwich".
#' @seealso \code{\link{adjust_loglik}} to adjust a user-supplied
#'   loglikelihood.
#' @export
coef.chandwich <- function(object, complete = TRUE, ...) {
  if (!inherits(object, "chandwich")) {
    stop("use only with \"chandwich\" objects")
  }
  if (complete) {
    cf <- attributes(small)$res_MLE
  } else {
    cf <- attributes(small)$MLE
  }
  return(cf)
}
