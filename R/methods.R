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
#'
#'   The default is \code{complete = FALSE}, which is in sync with
#'   \code{\link{vcov.chandwich}}.
#' @param ... Additional optional arguments. At present no optional
#'   arguments are used.
#' @details The full vector of estimates is taken from
#' \code{attributes(object)$res_MLE} and the reduced vector from
#' \code{attributes(object)$MLE}.
#' @return A numeric vector of estimated parameters, which will be
#'  named if names were provided in the call to
#'  \code{\link{adjust_loglik}}.
#' @seealso \code{\link{vcov.chandwich}}: \code{vcov} method for
#'   class "chandwich".
#' @seealso \code{\link{summary.chandwich}}: \code{summary} method for
#'   class "chandwich".
#' @seealso \code{\link{adjust_loglik}} to adjust a user-supplied
#'   loglikelihood.
#' @export
coef.chandwich <- function(object, complete = FALSE, ...) {
  if (!inherits(object, "chandwich")) {
    stop("use only with \"chandwich\" objects")
  }
  if (complete) {
    cf <- attributes(object)$res_MLE
  } else {
    cf <- attributes(object)$MLE
  }
  return(cf)
}

# ============================== vcov.chandwich ==============================

#' Calculate the variance-covariance model for an object of class "chandwich"
#'
#' \code{vcov} method for class "chandwich".
#'
#' @param object an object of class "chandwich", a result of a call to
#'   \code{\link{adjust_loglik}}.
#' @param complete A logical scalar.  If \code{complete = TRUE} then the
#'   full variance-covariance matrix is returned, including any fixed
#'   using \code{fixed_pars} and \code{fixed_at} in the call to
#'   \code{\link{adjust_loglik}}.  In this instance the corresponding
#'   elements of the matrix are \code{NA}. Otherwise, the variance-covariance
#'   matrix of only the free parameters is returned.
#'
#'   The default is \code{complete = FALSE}, which is in sync with
#'   \code{\link{coef.chandwich}}.
#' @param adjusted A logical scalar.  If \code{adjusted = TRUE} then the
#'   variance-covariance matrix is estimated using a sandwich estimator.
#'   Otherwise, the inverse of the observed information at the MLE is used.
#' @param ... Additional optional arguments. At present no optional
#'   arguments are used.
#' @details The variance-covariance matrix is based on
#' \code{attributes(object)$adjVC} for \code{adjusted = TRUE} and
#' \code{attributes(object)$VC} for \code{adjusted = FALSE}.
#' These return the estimate variance-covariance matrix of only the
#' free parameters.
#' @return A numeric matrix.  The dimensions will be named if names were
#'   provided in the call to \code{\link{adjust_loglik}}.
#' @seealso \code{\link{coef.chandwich}}: \code{coef} method for
#'   class "chandwich".
#' @seealso \code{\link{summary.chandwich}}: \code{summary} method for
#'   class "chandwich".
#' @seealso \code{\link{adjust_loglik}} to adjust a user-supplied
#'   loglikelihood.
#' @export
vcov.chandwich <- function(object, complete = FALSE, adjusted = TRUE, ...) {
  if (!inherits(object, "chandwich")) {
    stop("use only with \"chandwich\" objects")
  }
  if (adjusted) {
    vc <- attributes(object)$adjVC
  } else {
    vc <- attributes(object)$VC
  }
  if (complete) {
    np <- attributes(object)$p_full
    pfree <- attributes(object)$free_pars
    dummy <- matrix(NA, np, np)
    dummy[pfree, pfree] <- vc
    vc <- dummy
    full_par_names <- attributes(object)$full_par_names
    dimnames(vc) <- list(full_par_names, full_par_names)
  }
  return(vc)
}

# ============================= logLik.chandwich ==============================

#' Extract log-likelihood for objects of class "chandwich"
#'
#' \code{logLik} method for class "chandwich".
#'
#' @param object an object of class "chandwich", a result of a call to
#'   \code{\link{adjust_loglik}}.
#' @param ... Additional optional arguments. At present no optional
#'   arguments are used.
#' @details The value of the maximised (independence) loglikelihood
#' is extracted from \code{attr(object, "max_loglik")}.  It is also equal
#' to \code{object(attr(object, "MLE"))}.
#' @return An object of class \code{"logLik"}: a numeric scalar with
#' value equal to the maximised (independence) loglikelihood, that is,
#' the value of the independence loglikelihood at the MLE, and the
#' attribute \code{"df"}, which gives the number of free parameters estimated.
#' @seealso \code{\link{coef.chandwich}}: \code{coef} method for
#'   class "chandwich".
#' @seealso \code{\link{vcov.chandwich}}: \code{vcov} method for
#'   class "chandwich".
#' @seealso \code{\link{summary.chandwich}}: \code{summary} method for
#'   class "chandwich".
#' @seealso \code{\link{adjust_loglik}} to adjust a user-supplied
#'   loglikelihood.
#' @export
logLik.chandwich <- function(object, ...) {
  if (!inherits(object, "chandwich")) {
    stop("use only with \"chandwich\" objects")
  }
  val <- attr(object, "max_loglik")
  attr(val, "df") <- attr(object, "p_current")
  class(val) <- "logLik"
  return(val)
}
