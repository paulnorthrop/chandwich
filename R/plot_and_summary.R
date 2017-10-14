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
      y <- cbind(y, x(x_vals, type = "none"))
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

# ============================== plot.confint =================================

#' Plot diagnostics a confint object
#'
#' \code{plot} method for class "confint".
#'
#' @param x an object of class "confint", a result of a call to
#'   \code{\link{conf_intervals}}.
#' @param y Not used.
#' @param ... Additional arguments passed on to ...
#' @export
plot.confint <- function(x, y, ..., add_lines = TRUE) {
  if (!inherits(x, "confint")) {
    stop("use only with \"confint\" objects")
  }
  my_matplot <- function(x, y, ..., type) {
    matplot(x, y, type = "l", ...)
  }
  my_matplot(x$parameter_vals, x$prof_loglik_vals, ...)
  if (add_lines) {
    abline(h = x$cutoff)
    abline(v = x$prof_CI[, 1])
    abline(v = x$prof_CI[, 2])
  }
  return(invisible())
}

# check for equality of max_logliks?

# ============================== plot.confreg =================================

#' Plot diagnostics a confreg object
#'
#' \code{plot} method for class "confreg".
#'
#' @param x an object of class "confreg", a result of a call to
#'   \code{\link{conf_region}}.
#' @param y Not used.
#' @param ... Additional arguments passed to \code{\link[graphics]{contour}}.
#' @export
plot.confreg <- function(x, y = NULL, y2 = NULL, y3 = NULL, conf = 95, ...) {
  if (!inherits(x, "confreg")) {
    stop("use only with \"confreg\" objects")
  }
  x_range <- range(x$grid1, y$grid1, y2$grid1, y3$grid1, finite = TRUE)
  y_range <- range(x$grid2, y$grid2, y2$grid2, y3$grid2, finite = TRUE)
  # User-supplied arguments for contour.
  user_args <- list(...)
  # If labels is not supplied then set it to confidence level
  if (is.null(user_args$labels)) {
    user_args$labels <- conf
  }
  # If drawlabels is not supplied then make it FALSE unless conf has length > 1
  if (is.null(user_args$drawlabels)) {
    user_args$drawlabels <- length(conf) > 1
  }
  if (is.null(user_args$xlim)) {
    user_args$xlim <- x_range
  }
  if (is.null(user_args$ylim)) {
    user_args$ylim <- y_range
  }
  if (is.null(user_args$col)) {
    my_col <- 1:4
  } else {
    my_col <- rep(user_args$col, 4)
  }
  # Create plot using x
  max_loglik <- attr(x$object, "max_loglik")
  cutoff <- max_loglik - qchisq(conf  / 100, 2) / 2
  user_args$col <- my_col[1]
  for_contour <- c(list(x = x$grid1, y = x$grid2, z = x$prof_loglik,
                        levels = cutoff), user_args)
  do.call(graphics::contour, for_contour)
  # Add to plot using y
  if (!is.null(y)) {
    max_loglik <- attr(y$object, "max_loglik")
    cutoff <- max_loglik - qchisq(conf  / 100, 2) / 2
    user_args$col <- my_col[2]
    for_contour <- c(list(x = y$grid1, y = y$grid2, z = y$prof_loglik,
                          levels = cutoff, add = TRUE), user_args)
    do.call(graphics::contour, for_contour)
  }
  # Add to plot using y2
  if (!is.null(y2)) {
    max_loglik <- attr(y2$object, "max_loglik")
    cutoff <- max_loglik - qchisq(conf  / 100, 2) / 2
    user_args$col <- my_col[3]
    for_contour <- c(list(x = y2$grid1, y = y2$grid2, z = y2$prof_loglik,
                          levels = cutoff, add = TRUE), user_args)
    do.call(graphics::contour, for_contour)
  }
  # Add to plot using y3
  if (!is.null(y2)) {
    max_loglik <- attr(y3$object, "max_loglik")
    cutoff <- max_loglik - qchisq(conf  / 100, 2) / 2
    user_args$col <- my_col[4]
    for_contour <- c(list(x = y3$grid1, y = y3$grid2, z = y3$prof_loglik,
                          levels = cutoff, add = TRUE), user_args)
    do.call(graphics::contour, for_contour)
  }
  return(invisible())
}
