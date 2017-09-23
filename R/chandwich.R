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

#' Rat tumor data
#'
#' Tumor incidence in 71 groups of rate from Tarone (1982).
#' The matrix \code{rat} has 71 rows and 2 columns.
#' Each row relates to a different group of rats.
#' The first column (\code{y}) contains the number of rats with tumors.
#' The second column (\code{n}) contains the total number of rats.
#'
#' @format A matrix with 71 rows and 2 columns.
#' @source Table 5.1 of Gelman, A., Carlin, J. B., Stern, H. S. Dunson, D. B.,
#'  Vehtari, A. and Rubin, D. B. (2013) \emph{Bayesian Data Analysis},
#'  Chapman & Hall / CRC.
#'   \url{http://www.stat.columbia.edu/~gelman/book/data/rats.asc}
#' @references Tarone, R. E. (1982) The use of historical information in
#'   testing for a trend in proportions. \emph{Biometrics}, \strong{38},
#'   215-220. \url{https://doi.org/10.2307/2530304}
"rats"
