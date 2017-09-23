#' chandwich: Chandler-Bate Sandwich Loglikelihood Adjustment
#'
#' Performs adjustments of an independence loglikelihood using
#' a robust sandwich estimator of the parameter covariance matrix, based on
#' the methodology in
#' \href{http://dx.doi.org/10.1093/biomet/asm015}{Chandler and Beta (2007)}.
#' This can be used for cluster correlated data when interest lies in the
#' parameters of the marginal distributions.
#'
#' @details Add details
#'
#'   See \code{vignette("chandwich-vignette", package = "chandwich")} for an
#'   overview of the package.
#' @references Chandler, R. E. and Bate, S. (2007). Inference for clustered
#'   data using the independence loglikelihood. \emph{Biometrika},
#'   \strong{94}(1), 167-183. \url{http://dx.doi.org/10.1093/biomet/asm015}
#' @docType package
#' @name chandwich
#' @import methods
NULL

#' Oxford and Worthing annual maximum temperatures
#'
#' Annual maximum temperatures at Oxford and Worthing (England), for the
#' period 1901 to 1980.
#'
#' @format A dataframe with 80 rows and 2 columns, named Oxford and Worthing.
#' @source Tabony, R. C. (1983) Extreme value analysis in meteorology.
#'  \emph{The Meteorological Magazine}, \strong{112}, 77-98.
#' @references Chandler, R. E. and Bate, S. (2007). Inference for clustered
#'   data using the independence loglikelihood. \emph{Biometrika},
#'   \strong{94}(1), 167-183. \url{http://dx.doi.org/10.1093/biomet/asm015}
"owtemps"

