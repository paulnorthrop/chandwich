#' The Generalised Extreme Value Density Function
#'
#' Density function of the generalised extreme value (GEV) distribution.
#'
#' @param x Numeric vectors of quantiles.
#' @param loc,scale,shape Numeric vectors.
#'   Location, scale and shape parameters.
#'   All elements of \code{scale} must be positive.
#' @param log A logical scalar.  If TRUE the log-density is returned.
#' @details The distribution function of a GEV distribution with parameters
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
#' @return A numeric vector of value(s) of the density of the GEV distribution.
#'   The length of the result is determined by the maximum of the lengths
#'   of the numerical arguments for the other functions.
#'   The numerical arguments are recycled to the length of the result.
#' @references Jenkinson, A. F. (1955) The frequency distribution of the
#'   annual maximum (or minimum) of meteorological elements.
#'   \emph{Quart. J. R. Met. Soc.}, \strong{81}, 158-171.
#'   Chapter 3: \url{http://dx.doi.org/10.1002/qj.49708134804}
#' @references Coles, S. G. (2001) \emph{An Introduction to Statistical
#'   Modeling of Extreme Values}, Springer-Verlag, London.
#'   \url{http://dx.doi.org/10.1007/978-1-4471-3675-0_3}
#' @examples
#' gev_dens(-1:4, 1, 0.5, 0.8)
#' gev_dens(1:6, 1, 0.5, -0.2, log = TRUE)
#' gev_dens(1, shape = c(-0.2, 0.4))
#' @export
gev_dens <- function(x, loc = 0, scale = 1, shape = 0, log = FALSE) {
  if (any(scale <= 0)) {
    stop("invalid scale: scale must be positive.")
  }
  max_len <- max(length(x), length(loc), length(scale), length(shape))
  x <- rep_len(x, max_len)
  loc <- rep_len(loc, max_len)
  scale <- rep_len(scale, max_len)
  shape <- rep_len(shape, max_len)
  x <- (x - loc) / scale
  xx <- 1 + shape * x
  d <- ifelse(xx < 0, 0,
         ifelse(xx == 0 & shape == -1, 1,
           ifelse(xx == 0 & shape < -1, Inf,
             ifelse(abs(shape) > 1e-6,
               xx ^ (-(1 + 1 / shape)) * exp(-xx ^ (-1/ shape)),
      exp(-x + shape * x * (x - 2) / 2 - exp(-x + shape * x ^ 2 / 2))))))
  d <- d / scale
  if (log) {
    d <- log(d)
  }
  return(d)
}
