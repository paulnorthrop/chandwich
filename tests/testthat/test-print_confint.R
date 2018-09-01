context("print_confint")

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
conf_n <- conf_intervals(rat_res, type = "none")

check_same <- utils::capture.output(print(conf_v))
check_conf_v <- utils::capture.output(conf_v)
test_that("binom: print Ok for type = vertical", {
  testthat::expect_identical(check_same, check_conf_v)
})
check_same <- utils::capture.output(print(conf_n))
check_conf_n <- utils::capture.output(conf_n)
test_that("binom: print OK for type = none", {
  testthat::expect_identical(check_same, check_conf_n)
})

# --------- Misspecified Poisson model for negative binomial data ----------

# ... following Section 5.1 of the "Object-Oriented Computation of Sandwich
# Estimators" vignette of the sandwich package
# https://cran.r-project.org/web/packages/sandwich/vignettes/sandwich-OOP.pdf

# Simulate data
set.seed(123)
x <- stats::rnorm(250)
y <- stats::rnbinom(250, mu = exp(1 + x), size = 1)
# Fit misspecified Poisson model
fm_pois <- stats::glm(y ~ x + I(x^2), family = poisson)

# Contributions to the independence loglikelihood
pois_glm_loglik <- function(pars, y, x) {
  log_mu <- pars[1] + pars[2] * x + pars[3] * x ^ 2
  return(stats::dpois(y, lambda = exp(log_mu), log = TRUE))
}
pars <- c("alpha", "beta", "gamma")
# Linear model (gamma fixed at 0)
pois_lin <- adjust_loglik(pois_glm_loglik, y = y, x = x, par_names = pars,
                          fixed_pars = "gamma")

conf_1 <- conf_intervals(pois_lin, which_par = "alpha")
conf_2 <- conf_intervals(pois_lin, which_par = 1, type = "none")

check_same <- utils::capture.output(print(conf_1))
check_conf_1 <- utils::capture.output(conf_1)
test_that("pois: print OK for type = vertical", {
  testthat::expect_identical(check_same, check_conf_1)
})
check_same <- utils::capture.output(print(conf_2))
check_conf_2 <- utils::capture.output(conf_2)
test_that("pois: print OK for type = none", {
  testthat::expect_identical(check_same, check_conf_2)
})
