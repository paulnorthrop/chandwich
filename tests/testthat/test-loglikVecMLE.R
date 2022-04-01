#context("loglikVecMLE")

# Check that the sum of the loglikelihood contributions at the MLE is
# equal to the value of max_loglik

# ------------------------- Binomial model, rats data ----------------------

# Contributions to the independence loglikelihood
binom_loglik <- function(prob, data) {
  if (prob < 0 || prob > 1) {
    return(-Inf)
  }
  return(stats::dbinom(data[, "y"], data[, "n"], prob, log = TRUE))
}

rat_res <- adjust_loglik(loglik = binom_loglik, data = rats, par_names = "p")

max_loglik <- attr(rat_res, "max_loglik")
max_loglik2 <- sum(attr(rat_res, "loglikVecMLE"))

test_that("rats: max_loglik and sum(loglikVecMLE) agree", {
  testthat::expect_identical(max_loglik, max_loglik2)
})

# ----------------------- GEV, Oxford and Worthing data -----------------------

# GEV independence loglikelihood for the Oxford-Worthing annual maximum
# temperature dataset owtemps

gev_loglik <- function(pars, data) {
  o_pars <- pars[c(1, 3, 5)] + pars[c(2, 4, 6)]
  w_pars <- pars[c(1, 3, 5)] - pars[c(2, 4, 6)]
  if (isTRUE(o_pars[2] <= 0 | w_pars[2] <= 0)) return(-Inf)
  o_data <- data[, "Oxford"]
  w_data <- data[, "Worthing"]
  check <- 1 + o_pars[3] * (o_data - o_pars[1]) / o_pars[2]
  if (isTRUE(any(check <= 0))) return(-Inf)
  check <- 1 + w_pars[3] * (w_data - w_pars[1]) / w_pars[2]
  if (isTRUE(any(check <= 0))) return(-Inf)
  o_loglik <- log_gev(o_data, o_pars[1], o_pars[2], o_pars[3])
  w_loglik <- log_gev(w_data, w_pars[1], w_pars[2], w_pars[3])
  return(o_loglik + w_loglik)
}

# Initial estimates (method of moments for the Gumbel case)
sigma <- as.numeric(sqrt(6 * diag(stats::var(owtemps))) / pi)
mu <- as.numeric(colMeans(owtemps) - 0.57722 * sigma)
init <- c(mean(mu), -diff(mu) / 2, mean(sigma), -diff(sigma) / 2, 0, 0)
par_names <-  c("mu0", "mu1", "sigma0", "sigma1", "xi0", "xi1")

# Full model
large <- adjust_loglik(gev_loglik, data = owtemps, init = init,
                       par_names = par_names)

max_loglik <- attr(large, "max_loglik")
max_loglik2 <- sum(attr(large, "loglikVecMLE"))

test_that("GEV large: max_loglik and sum(loglikVecMLE) agree", {
  testthat::expect_identical(max_loglik, max_loglik2)
})

# Restricted model, using larger to start and character name
medium <- adjust_loglik(larger = large, fixed_pars = "xi1")

max_loglik <- attr(medium, "max_loglik")
max_loglik2 <- sum(attr(medium, "loglikVecMLE"))

test_that("GEV medium: max_loglik and sum(loglikVecMLE) agree", {
  testthat::expect_identical(max_loglik, max_loglik2)
})

small <- adjust_loglik(larger = medium, fixed_pars = c("sigma1", "xi1"))

max_loglik <- attr(small, "max_loglik")
max_loglik2 <- sum(attr(small, "loglikVecMLE"))

test_that("GEV small: max_loglik and sum(loglikVecMLE) agree", {
  testthat::expect_identical(max_loglik, max_loglik2)
})

