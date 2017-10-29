# ============================== plot.chandwich ===============================

#' Plot diagnostics a chandwich object
#'
#' \code{plot} method for class "chandwich".  Only applicable to an object
#' \code{x} for which \code{attr(x, "p_current") = 1}, i.e. a model with
#' one free parameter.
#'
#' @param x an object of class "chandwich", a result of a call to
#'   \code{\link{adjust_loglik}}.
#' @param y Not used.
#' @param type An integer vector, a subset of the numbers \code{1:4}.
#'   Indicates which loglikelihoods to plot: \code{1} for \code{"vertical"}
#'   adjustment; \code{2} for \code{"cholesky"} (horizontal adjustment);
#'   \code{3} for \code{"spectral"} (horizontal adjustment); \code{4}
#'   for no adjustment, i.e. based on the independence loglikelihood.
#' @param legend A logical scalar or a character vector.  If this is
#'   supplied then a legend is added to the plot.  If \code{legend} is a
#'   character vector then it is used as the argument \code{legend}
#'   to \code{\link[graphics]{legend}}.  Otherwise, i.e. if
#'   \code{legend = TRUE} then the argument \code{type} is used.
#' @param legend_pos The position of the legend (if required) specified using
#'   the argument \code{x} in \code{\link[graphics]{legend}}.
#' @param ... Additional arguments passed to \code{\link[graphics]{matplot}}
#'   or \code{\link[graphics]{legend}}.  The arguments \code{col}, \code{lty}
#'   and \code{lwd} will (in a consistent way) by both
#'   \code{\link[graphics]{matplot}} and \code{\link[graphics]{legend}}.
#'
#'   If the argument \code{xlim} to \code{\link[graphics]{matplot}} is not
#'   supplied then the MLE minus (for \code{lower}) or plus (for \code{upper})
#'   standard errors is used.  If \code{type} does not include 4 then adjusted
#'   standard errors are used.  Otherwise, the larger of the adjust and
#'   unadjusted standard errors are used.
#' @return Nothing is returned.
#' @seealso \code{\link{conf_region}} and \code{\link{plot.confreg}} to
#'   plot a confidence region for a pair of parameters.
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
#' # Vertically adjusted loglikelihood only
#' plot(rat_res)
#' # Three adjusted loglikelihoods and the independence loglikelihood
#' plot(rat_res, type = 1:4)
#' # Plot over (0,1) and reposition the legend
#' plot(rat_res, type = 1:4, xlim = c(0, 1), legend_pos = "bottom")
#' @seealso \code{\link{adjust_loglik}} to adjust a user-supplied
#'   loglikelhood function.
#' @seealso \code{\link{summary.chandwich}} for maximum likelihood estimates
#'   and unadjusted and adjusted standard errors.
#' @export
plot.chandwich <- function(x, y, type = 1, legend = length(type) > 1,
                           legend_pos = "topleft", ...) {
  if (!inherits(x, "chandwich")) {
    stop("use only with \"chandwich\" objects")
  }
  n_pars <- attr(x, "p_current")
  if (n_pars != 1) {
    stop("x must have one free parameter, i'e. attr(x, ''p_current'') = 1")
  }
  # User-supplied arguments
  user_args <- list(...)
  # Always plot lines only
  user_args$type = "l"
  l_cond <- names(user_args) %in% methods::formalArgs(graphics::legend)
  lines_cond <- names(user_args) %in% c("col", "lty", "lwd")
  legend_args <- user_args[l_cond]
  user_args <- user_args[!l_cond | lines_cond]
  # If xlab or ylab are not supplied then use names(x$which_pars), if present
  if (is.null(user_args$xlab)) {
    if (!is.null(names(attr(x, "free_pars")))) {
      xlabel <- names(attr(x, "free_pars"))
    } else {
      xlabel <- ""
    }
    user_args$xlab <- parse(text = xlabel)
  }
  if (is.null(user_args$ylab)) {
    if (attr(x, "p_current") == 1) {
      user_args$ylab <- "loglikelihood"
    } else {
      user_args$ylab <- "profile loglikelihood"
    }
  }
  if (is.null(user_args$xlim)) {
    if (4 %in% type) {
      lower <- attr(x, "MLE") - 3 * max(attr(x, "adjSE"), attr(x, "SE"))
      upper <- attr(x, "MLE") + 3 * max(attr(x, "adjSE"), attr(x, "SE"))
    } else {
      lower <- attr(x, "MLE") - 3 * attr(x, "adjSE")
      upper <- attr(x, "MLE") + 3 * attr(x, "adjSE")
    }
    user_args$xlim <- c(lower, upper)
  }
  if (is.null(user_args$col)) {
    user_args$col <- rep(1, 4)
  }
  legend_args$col <- user_args$col
  if (is.null(user_args$lty)) {
    user_args$lty <- 1:4
  }
  legend_args$lty <- user_args$lty
  if (is.null(user_args$lwd)) {
    user_args$lwd <- rep(1, 4)
  }
  legend_args$lwd <- user_args$lwd
  # Create values for the plot
  x_vals <- seq(user_args$xlim[1], user_args$xlim[2], len = 100)
  y_vals <- NULL
  if (any(type == 1)) {
    y_vals <- cbind(y_vals, x(x_vals, type = "vertical"))
  }
  if (any(type == 2)) {
    y_vals <- cbind(y_vals, x(x_vals, type = "cholesky"))
  }
  if (any(type == 3)) {
    y_vals <- cbind(y_vals, x(x_vals, type = "spectral"))
  }
  if (any(type == 4)) {
    y_vals <- cbind(y_vals, x(x_vals, type = "none"))
  }
  for_matplot <- c(list(x = x_vals, y = y_vals), user_args)
  do.call(graphics::matplot, for_matplot)
  # Add a legend?
  if (legend | is.character(legend)) {
    types <- c("vertical", "cholesky", "spectral", "none")
    legend_args$x <- legend_pos
    if (legend) {
      legend_args$legend <- types[type]
    } else {
      legend_args$legend <- legend
    }
    if (is.null(legend_args$title)) {
      legend_args$title <- "adjustment"
    }
    do.call(graphics::legend, legend_args)
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
#' @seealso \code{\link{adjust_loglik}} to adjust a user-supplied
#'   loglikelhood function.
#' @seealso \code{\link{plot.chandwich}} for plots of one-dimensional adjusted
#'   loglikelihoods.
#' @section Examples:
#' See the examples in \code{\link{adjust_loglik}}.
#' @export
summary.chandwich <- function(object, digits = max(3, getOption("digits")-3),
                              ...) {
  if (!inherits(object, "chandwich")) {
    stop("use only with \"chandwich\" objects")
  }
  mle <- signif(attr(object, "MLE"), digits = digits)
  SE <- signif(attr(object, "SE"), digits = digits)
  adjSE <- signif(attr(object, "adjSE"), digits = digits)
  res <- cbind(MLE = mle, SE = SE, `adj. SE` = adjSE)
  if (is.null(attr(object, "par_names"))) {
    rownames(res) <- 1:length(mle)
  }
  return(res)
}

# ============================== plot.confreg =================================

#' Plot diagnostics a confreg object
#'
#' \code{plot} method for class "confreg".
#' Plots confidence regions for pairs of parameters using the profile
#' loglikelihood values calculated by \code{\link{conf_region}}.
#' Up to 4 different types of loglikelihood (see the argument \code{type}
#' to the function returned by \code{\link{adjust_loglik}})
#' may be superimposed on the same plot.
#'
#' @param x,y,y2,y3 objects of class "confreg", results of calls to
#'   \code{\link{conf_region}} for a common model and a common value of
#'   \code{which_pars}.  Contours are plotted for each object.
#' @param conf A numeric vector of confidence levels, i.e. numbers in
#'   (0, 100).  A confidence region contour is plotted for each value in
#'   \code{conf}.
#' @param legend A logical scalar or a character vector.  If this is
#'   supplied then a legend is added to the plot.  If \code{legend} is a
#'   character vector then it is used as the argument \code{legend}
#'   to \code{\link[graphics]{legend}}.  Otherwise, i.e. if
#'   \code{legend = TRUE} then the component \code{type} of the input
#'   object(s) \code{x, y, y2, y3} are used.
#' @param legend_pos The position of the legend (if required) specified using
#'   the argument \code{x} in \code{\link[graphics]{legend}}.
#' @param ... Additional arguments passed to \code{\link[graphics]{contour}}
#'  or \code{\link[graphics]{legend}}.  The arguments \code{col}, \code{lty}
#'  and \code{lwd} will (in a consistent way) by both
#'  \code{\link[graphics]{contour}} and \code{\link[graphics]{legend}}.
#' @return Nothing is returned.
#' @section Examples:
#' See the examples in \code{\link{conf_region}}.
#' @export
plot.confreg <- function(x, y = NULL, y2 = NULL, y3 = NULL, conf = 95,
                         legend = any(c(!is.null(y), !is.null(y2),
                                        !is.null(y3))),
                         legend_pos = "topleft", ...) {
  if (!inherits(x, "confreg")) {
    stop("x must be a \"confreg\" object")
  }
  check_name <- x$name
  # Check that all supplied objects have profiled the same parameters
  if (!is.null(y)) {
    if (!inherits(y, "confreg")) {
      stop("y must be a \"confreg\" object")
    }
    if (!identical(x$which_pars, y$which_pars)) {
      stop("y$which_pars is not identical to x$which_pars")
    }
    if (!identical(check_name, y$name)) {
      stop("y is not derived from the same model as x")
    }
  }
  if (!is.null(y2)) {
    if (!inherits(y2, "confreg")) {
      stop("y2 must be a \"confreg\" object")
    }
    if (!identical(x$which_pars, y2$which_pars)) {
      stop("y2$which_pars is not identical to x$which_pars")
    }
    if (!identical(check_name, y2$name)) {
      stop("y2 is not derived from the same model as x")
    }
  }
  if (!is.null(y3)) {
    if (!inherits(y3, "confreg")) {
      stop("y3 must be a \"confreg\" object")
    }
    if (!identical(x$which_pars, y3$which_pars)) {
      stop("y3$which_pars is not identical to x$which_pars")
    }
    if (!identical(check_name, y3$name)) {
      stop("y3 is not derived from the same model as x")
    }
  }
  # User-supplied arguments for contour
  user_args <- list(...)
  l_cond <- names(user_args) %in% methods::formalArgs(graphics::legend)
  lines_cond <- names(user_args) %in% c("col", "lty", "lwd")
  legend_args <- user_args[l_cond]
  user_args <- user_args[!l_cond | lines_cond]
  # If xlab or ylab are not supplied then use names(x$which_pars), if present
  if (is.null(user_args$xlab)) {
    user_args$xlab <- parse(text = names(x$which_pars)[1])
  }
  if (is.null(user_args$ylab)) {
    user_args$ylab <- parse(text = names(x$which_pars)[2])
  }
  # If labels is not supplied then set it to confidence level
  if (is.null(user_args$labels)) {
    user_args$labels <- conf
  }
  # If drawlabels is not supplied then make it FALSE unless conf has length > 1
  if (is.null(user_args$drawlabels)) {
    user_args$drawlabels <- length(conf) > 1
  }
  if (is.null(user_args$xlim)) {
    x_range <- range(x$grid1, y$grid1, y2$grid1, y3$grid1, finite = TRUE)
    user_args$xlim <- x_range
  }
  if (is.null(user_args$ylim)) {
    y_range <- range(x$grid2, y$grid2, y2$grid2, y3$grid2, finite = TRUE)
    user_args$ylim <- y_range
  }
  if (is.null(user_args$col)) {
    my_col <- rep(1, 4)
  } else {
    my_col <- rep_len(user_args$col, 4)
  }
  legend_args$col <- my_col
  if (is.null(user_args$lty)) {
    my_lty <- 1:4
  } else {
    my_lty <- rep_len(user_args$lty, 4)
  }
  legend_args$lty <- my_lty
  if (is.null(user_args$lwd)) {
    my_lwd <- rep(1, 4)
  } else {
    my_lwd <- rep_len(user_args$lwd, 4)
  }
  legend_args$lwd <- my_lwd
  # Create plot using x
  max_loglik <- x$max_loglik
  cutoff <- max_loglik - stats::qchisq(conf  / 100, 2) / 2
  user_args$col <- my_col[1]
  user_args$lty <- my_lty[1]
  user_args$lwd <- my_lwd[1]
  for_contour <- c(list(x = x$grid1, y = x$grid2, z = x$prof_loglik,
                        levels = cutoff), user_args)
  do.call(graphics::contour, for_contour)
  types <- x$type
  # Add to plot using y
  if (!is.null(y)) {
    max_loglik <- y$max_loglik
    cutoff <- max_loglik - stats::qchisq(conf  / 100, 2) / 2
    user_args$col <- my_col[2]
    user_args$lty <- my_lty[2]
    user_args$lwd <- my_lwd[2]
    for_contour <- c(list(x = y$grid1, y = y$grid2, z = y$prof_loglik,
                          levels = cutoff, add = TRUE), user_args)
    do.call(graphics::contour, for_contour)
    types <- c(types, y$type)
  }
  # Add to plot using y2
  if (!is.null(y2)) {
    max_loglik <- y2$max_loglik
    cutoff <- max_loglik - stats::qchisq(conf  / 100, 2) / 2
    user_args$col <- my_col[3]
    user_args$lty <- my_lty[3]
    user_args$lwd <- my_lwd[3]
    for_contour <- c(list(x = y2$grid1, y = y2$grid2, z = y2$prof_loglik,
                          levels = cutoff, add = TRUE), user_args)
    do.call(graphics::contour, for_contour)
    types <- c(types, y2$type)
  }
  # Add to plot using y3
  if (!is.null(y3)) {
    max_loglik <- y3$max_loglik
    cutoff <- max_loglik - stats::qchisq(conf  / 100, 2) / 2
    user_args$col <- my_col[4]
    user_args$lty <- my_lty[4]
    user_args$lwd <- my_lwd[4]
    for_contour <- c(list(x = y3$grid1, y = y3$grid2, z = y3$prof_loglik,
                          levels = cutoff, add = TRUE), user_args)
    do.call(graphics::contour, for_contour)
    types <- c(types, y3$type)
  }
  # Add a legend?
  if (legend | is.character(legend)) {
    legend_args$x <- legend_pos
    if (legend) {
      legend_args$legend <- types
    } else {
      legend_args$legend <- legend
    }
    if (is.null(legend_args$title)) {
      legend_args$title <- "adjustment"
    }
    do.call(graphics::legend, legend_args)
  }
  return(invisible())
}

# ============================== plot.confint =================================

#' Plot diagnostics a confint object
#'
#' \code{plot} method for class "confint".
#' Plots the (profile) loglikelihood for a parameter using the values
#' calculated by \code{\link{conf_intervals}}.
#' Up to 4 different types of loglikelihood (see the argument \code{type}
#' to the function returned by \code{\link{adjust_loglik}})
#' may be superimposed on the same plot.
#' By default (\code{add_lines = TRUE}) the confidence limits calculated
#' in \code{\link{conf_intervals}} are indicated on the plot .
#'
#' @param x,y,y2,y3 objects of class "confint", results of calls to
#'   \code{\link{conf_intervals}} for a common model.  A (profile)
#'   loglikelihood will be plotted for each object.
#' @param which_par A scalar specifying the parameter for which the plot
#'   is produced.  Can be either a numeric vector, specifying index of the
#'   component of the \strong{full} parameter vector, or a character scalar
#'   parameter name.  The former must be in \code{x$which_pars}, the latter
#'   must be in \code{names(x$which_pars)}.
#' @param add_lines A logical scalar.  Whether or not to add horizontal lines
#'   to the plot to identify the confidence limits.
#' @param conf A numeric vector of values in (0, 100).  If
#'   \code{add_lines = TRUE} then a horizontal line is added for each
#'   value in \code{conf}.  If \code{conf} is not supplied then the
#'   value stored in \code{x$conf} is used.
#' @param legend A logical scalar or a character vector.  If this is
#'   supplied then a legend is added to the plot.  If \code{legend} is a
#'   character vector then it is used as the argument \code{legend}
#'   to \code{\link[graphics]{legend}}.  Otherwise, i.e. if
#'   \code{legend = TRUE} then the component \code{type} of the input
#'   object(s) \code{x, y, y2, y3} are used.
#' @param legend_pos The position of the legend (if required) specified using
#'   the argument \code{x} in \code{\link[graphics]{legend}}.
#' @param ... Additional arguments passed to \code{\link[graphics]{matplot}}
#'   or \code{\link[graphics]{legend}}.  The arguments \code{col}, \code{lty}
#'   and \code{lwd} will (in a consistent way) by both
#'   \code{\link[graphics]{matplot}} and \code{\link[graphics]{legend}}.
#' @return Nothing is returned.
#' @section Examples:
#' See the examples in \code{\link{conf_intervals}}.
#' @export
plot.confint <- function(x, y = NULL, y2 = NULL, y3 = NULL,
                         which_par = x$which_pars[1], conf = x$conf,
                         add_lines = TRUE,
                         legend = any(c(!is.null(y), !is.null(y2),
                                        !is.null(y3))),
                         legend_pos = "topleft", ...) {
  if (!inherits(x, "confint")) {
    stop("x must be a \"confint\" object")
  }
  if (length(which_par) != 1) {
    stop("which_par must have length one")
  }
  # Where is which_par positioned in x$which_pars
  if (is.numeric(which_par)) {
    which_index <- match(which_par, x$which_pars)
    if (!is.null(names(x$which_pars))) {
      par_names <- names(x$which_pars)
      xlabel <- par_names[which_index]
    } else {
      xlabel <- ""
    }
  } else if (is.character(which_par)) {
    if (is.null(names(x$which_pars))) {
      stop("which_par can be character only if names(x$which_pars) isn't NULL")
    }
    par_names <- names(x$which_pars)
    xlabel <- par_names[which_index]
    which_index <- match(which_par, par_names)
    which_par <- x$which_pars[which_index]
  } else {
    stop("which_par must be numeric or character")
  }
  # Function to extract the desired values from the input objects
  extract_values <- function(object) {
    if (which_par %in% object$which_pars) {
      which_column <- match(which_par, object$which_pars)
      x_vals <- object$parameter_vals[, which_column]
      y_vals <- object$prof_loglik_vals[, which_column]
    } else {
      stop("Numeric which_par must be in object$which_pars")
    }
    return(list(x_vals = x_vals, y_vals = y_vals))
  }
  # Pick the correct columns in x$parameter_vals and x$prof_loglik_vals
  temp <- extract_values(x)
  x_vals <- temp$x_vals
  y_vals <- temp$y_vals
  types <- x$type
  check_name <- attr(x, "name")
  if (!is.null(y)) {
    if (!inherits(y, "confint")) {
      stop("y must be a \"confint\" object")
    }
    temp <- extract_values(y)
    x_vals <- cbind(x_vals, temp$x_vals)
    y_vals <- cbind(y_vals, temp$y_vals)
    types <- c(types, y$type)
  }
  if (!is.null(y2)) {
    if (!inherits(y2, "confint")) {
      stop("y2 must be a \"confint\" object")
    }
    temp <- extract_values(y2)
    x_vals <- cbind(x_vals, temp$x_vals)
    y_vals <- cbind(y_vals, temp$y_vals)
    types <- c(types, y2$type)
  }
  if (!is.null(y3)) {
    if (!inherits(y3, "confint")) {
      stop("y3 must be a \"confint\" object")
    }
    temp <- extract_values(y3)
    x_vals <- cbind(x_vals, temp$x_vals)
    y_vals <- cbind(y_vals, temp$y_vals)
    types <- c(types, y3$type)
  }
  # User-supplied arguments
  user_args <- list(...)
  # Always plot lines only
  user_args$type = "l"
  l_cond <- names(user_args) %in% methods::formalArgs(graphics::legend)
  lines_cond <- names(user_args) %in% c("col", "lty", "lwd")
  legend_args <- user_args[l_cond]
  user_args <- user_args[!l_cond | lines_cond]
  # If xlab or ylab are not supplied then use names(x$which_pars), if present
  if (is.null(user_args$xlab)) {
    user_args$xlab <- parse(text = xlabel)
  }
  if (is.null(user_args$ylab)) {
    if (x$p_current == 1) {
      user_args$ylab <- "loglikelihood"
    } else {
      user_args$ylab <- "profile loglikelihood"
    }
  }
  if (is.null(user_args$xlim)) {
    x_range <- range(x_vals, finite = TRUE)
    user_args$xlim <- x_range
  }
  if (is.null(user_args$ylim)) {
    y_range <- range(y_vals, finite = TRUE)
    user_args$ylim <- y_range
  }
  if (is.null(user_args$col)) {
    user_args$col <- rep(1, 4)
  }
  legend_args$col <- user_args$col
  if (is.null(user_args$lty)) {
    user_args$lty <- 1:4
  }
  legend_args$lty <- user_args$lty
  if (is.null(user_args$lwd)) {
    user_args$lwd <- rep(1, 4)
  }
  legend_args$lwd <- user_args$lwd
  for_matplot <- c(list(x = x_vals, y = y_vals), user_args)
  do.call(graphics::matplot, for_matplot)
  cutoff <- x$max_loglik - stats::qchisq(conf  / 100, 1) / 2
  if (add_lines) {
    graphics::abline(h = cutoff)
  }
  # Add a legend?
  if (legend | is.character(legend)) {
    legend_args$x <- legend_pos
    if (legend) {
      legend_args$legend <- types
    } else {
      legend_args$legend <- legend
    }
    if (is.null(legend_args$title)) {
      legend_args$title <- "adjustment"
    }
    do.call(graphics::legend, legend_args)
  }
  return(invisible())
}
