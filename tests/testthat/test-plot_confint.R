context("conf_intervals")

# -------------------------- GEV model, owtemps data -----------------------
# ------------ following Section 5.2 of Chandler and Bate (2007) -----------

gev_loglik <- function(pars, data) {
  o_pars <- pars[c(1, 3, 5)] + pars[c(2, 4, 6)]
  w_pars <- pars[c(1, 3, 5)] - pars[c(2, 4, 6)]
  if (o_pars[2] <= 0 | w_pars[2] <= 0) return(-Inf)
  o_data <- data[, "Oxford"]
  w_data <- data[, "Worthing"]
  check <- 1 + o_pars[3] * (o_data - o_pars[1]) / o_pars[2]
  if (any(check <= 0)) return(-Inf)
  check <- 1 + w_pars[3] * (w_data - w_pars[1]) / w_pars[2]
  if (any(check <= 0)) return(-Inf)
  o_loglik <- log_gev(o_data, o_pars[1], o_pars[2], o_pars[3])
  w_loglik <- log_gev(w_data, w_pars[1], w_pars[2], w_pars[3])
  return(o_loglik + w_loglik)
}

# Initial estimates (method of moments for the Gumbel case)
sigma <- as.numeric(sqrt(6 * diag(var(owtemps))) / pi)
mu <- as.numeric(colMeans(owtemps) - 0.57722 * sigma)
init <- c(mean(mu), -diff(mu) / 2, mean(sigma), -diff(sigma) / 2, 0, 0)

# Log-likelihood adjustment of the full model
par_names <- c("mu[0]", "mu[1]", "sigma[0]", "sigma[1]", "xi[0]", "xi[1]")
large <- adjust_loglik(gev_loglik, data = owtemps, init = init,
                       par_names = par_names)

# 95% likelihood-based confidence intervals, vertically adjusted
large_v <- conf_intervals(large, which_pars = "xi[1]")

check_NULL <- try(plot(large_v), silent = TRUE)
test_that("Default which_par gives no error", {
  testthat::expect_identical(check_NULL, NULL)
})
check_NULL <- try(plot(large_v, which_par = 6), silent = TRUE)
test_that("Default which_par = 6 gives no error", {
  testthat::expect_identical(check_NULL, NULL)
})
check_NULL <- try(plot(large_v, which_par = "xi[1]"), silent = TRUE)
test_that("Default which_par = xi[1], gives no error", {
  testthat::expect_identical(check_NULL, NULL)
})

check_error <- try(plot(large_v, which_par = 1), silent = TRUE)
test_that("Inappropriate which_par gives an error", {
  testthat::expect_identical(class(check_error), "try-error")
})

# ------------------------- Binomial model, rats data ----------------------

# Contributions to the independence loglikelihood
binom_loglik <- function(prob, data) {
  if (prob < 0 || prob > 1) {
    return(-Inf)
  }
  return(dbinom(data[, "y"], data[, "n"], prob, log = TRUE))
}
rat_res <- adjust_loglik(loglik = binom_loglik, data = rats, par_names = "p")

conf_v <- conf_intervals(rat_res)
conf_none <- conf_intervals(rat_res, type = "none")
conf_s <- conf_intervals(rat_res, type = "spectral")
conf_c <- conf_intervals(rat_res, type = "cholesky")

check_NULL <- try(plot(conf_v, conf_none, conf_s, conf_c), silent = TRUE)
test_that("Plotting 4 intervals gives no error", {
  testthat::expect_identical(check_NULL, NULL)
})

check_error <- try(plot(conf_v, conf_none, conf_s, conf_c, which_par = 3),
                   silent = TRUE)
test_that("Wrong parameter number gives an error", {
  testthat::expect_identical(class(check_error), "try-error")
})

