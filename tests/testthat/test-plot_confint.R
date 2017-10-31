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
