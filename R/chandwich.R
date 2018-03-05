#' chandwich: Chandler-Bate Sandwich Loglikelihood Adjustment
#'
#' Performs adjustments of an independence loglikelihood using
#' a robust sandwich estimator of the parameter covariance matrix, based on
#' the methodology in
#' \href{http://dx.doi.org/10.1093/biomet/asm015}{Chandler and Bate (2007)}.
#' This can be used for cluster correlated data when interest lies in the
#' parameters of the marginal distributions.
#' Functions for profiling the adjusted loglikelihoods are also provided, as
#' are functions for calculating and plotting confidence intervals, for single
#' model parameters, and confidence regions, for pairs of model parameters.
#'
#' @details
#' The main function in the chandwich package is \code{adjust_loglik}.  It
#' finds the maximum likelihood estimate (MLE) of model parameters based on
#' an independence loglikelihood in which cluster dependence in the data is
#' ignored.  The independence loglikelihood is adjusted in a way that ensures
#' that the Hessian of the adjusted loglikelihood coincides with a robust
#' sandwich estimate of the parameter covariance at the MLE.  Three
#' adjustments are available: one in which the independence loglikelihood
#' itself is scaled (vertical scaling) and two others where the scaling
#' is in the parameter vector (horizontal scaling).
#'
#' See Chandler and Bate (2007) for full details and
#' \code{vignette("chandwich-vignette", package = "chandwich")} for an
#' overview of the package.
#' @references Chandler, R. E. and Bate, S. (2007). Inference for clustered
#'   data using the independence loglikelihood. \emph{Biometrika},
#'   \strong{94}(1), 167-183. \url{http://dx.doi.org/10.1093/biomet/asm015}
#' @seealso \code{\link{adjust_loglik}} to adjust the independence
#'   loglikelihood.
#' @seealso \code{\link{compare_models}} to compare nested models using an
#'   adjusted loglikelihood ratio test.
#' @seealso \code{\link{conf_intervals}} to calculate confidence intervals
#'   for individual model parameters.
#' @seealso \code{\link{conf_region}} to calculate a confidence region
#'   for a pair of model parameters.
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
