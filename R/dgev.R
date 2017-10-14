#' The Generalised Extreme Value Log-Density Function
#'
#' Log-Density function of the generalised extreme value (GEV) distribution
#'
#' @param x Numeric vectors of quantiles.
#' @param loc,scale,shape Numeric scalars.
#'   Location, scale and shape parameters.
#'   \code{scale} must be positive.
#' @details \strong{It is assumed that \code{x}, \code{loc} = \eqn{\mu},
#'  \code{scale} = \eqn{\sigma} and \code{shape} = \eqn{\xi} are such that
#'  the GEV density is non-zero, i.e. that
#'  \eqn{1 + \xi (x - \mu) / \sigma > 0}. No check of this, or that
#'  \code{scale} > 0 is performed in this function.}
#'
#'  The distribution function of a GEV distribution with parameters
#'  \code{loc} = \eqn{\mu}, \code{scale} = \eqn{\sigma} (>0) and
#'  \code{shape} = \eqn{\xi} is
#'  \deqn{F(x) = exp { - [1 + \xi (x - \mu) / \sigma] ^ (-1/\xi)} }
#'  for \eqn{1 + \xi (x - \mu) / \sigma > 0}.  If \eqn{\xi = 0} the
#'  distribution function is defined as the limit as \eqn{\xi} tends to zero.
#'  The support of the distribution depends on \eqn{\xi}: it is
#'  \eqn{x <= \mu - \sigma / \xi} for \eqn{\xi < 0};
#'  \eqn{x >= \mu - \sigma / \xi} for \eqn{\xi > 0};
#'  and \eqn{x} is unbounded for \eqn{\xi = 0}.
#'  Note that if \eqn{\xi < -1} the GEV density function becomes infinite
#'  as \eqn{x} approaches \eqn{\mu -\sigma / \xi} from below.
#'
#'  See
#'  \url{https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution}
#'  for further information.
#'
#' @return A numeric vector of value(s) of the log-density of the GEV distribution.
#' @references Jenkinson, A. F. (1955) The frequency distribution of the
#'   annual maximum (or minimum) of meteorological elements.
#'   \emph{Quart. J. R. Met. Soc.}, \strong{81}, 158-171.
#'   Chapter 3: \url{http://dx.doi.org/10.1002/qj.49708134804}
#' @references Coles, S. G. (2001) \emph{An Introduction to Statistical
#'   Modeling of Extreme Values}, Springer-Verlag, London.
#'   \url{http://dx.doi.org/10.1007/978-1-4471-3675-0_3}
#' @examples
#' log_gev(1:4, 1, 0.5, 0.8)
#' log_gev(1:3, 1, 0.5, -0.2)
#' @export
log_gev <- function(x, loc = 0, scale = 1, shape = 0) {
  x <- (x - loc) / scale
  xx <- 1 + shape * x
  if (abs(shape) > 1e-6) {
    logd <- -(1 + 1 / shape) * log(xx) - xx ^ (-1 / shape)
  } else {
    logd <- -x + shape * x * (x - 2) / 2 - exp(-x + shape * x ^ 2 / 2)
  }
  return(logd - log(scale))
}
