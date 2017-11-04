context("conf_region and print")

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
pois_v <- conf_region(pois_lin)
pois_n <- conf_region(pois_lin, type = "none")
pois_c <- conf_region(pois_lin, type = "cholesky")
pois_s <- conf_region(pois_lin, type = "spectral")
check_NULL <- try(plot(pois_n, pois_v, pois_c, pois_s,
                       conf = c(50, 95),
                       col = 4:1, lwd = 2, lty = 1), silent = TRUE)

test_that("Consistent names gives no error", {
  testthat::expect_identical(check_NULL, NULL)
})

# Repeat for character which_pars
pois_v_char <- conf_region(pois_lin, which_pars = c("alpha", "beta"))
check_NULL <- try(plot(pois_n, pois_v_char, pois_c, pois_s), silent = TRUE)

test_that("Consistent names gives no error, which_pars is character", {
  testthat::expect_identical(check_NULL, NULL)
})

new_pois_lin <- adjust_loglik(pois_glm_loglik, y = y, x = x, par_names = pars,
                              fixed_pars = "gamma", name = "wrong_name")
new_pois_n <- conf_region(new_pois_lin, type = "none")
check_error <- try(plot(new_pois_n, pois_v, conf = c(50, 95)),
                   silent = TRUE)

test_that("Inconsistent names gives an error", {
  testthat::expect_identical(class(check_error), "try-error")
})

# Print

check_same <- try(print(pois_v), silent = TRUE)
test_that("Printing gives no error for type = vertical", {
  testthat::expect_identical(check_same, pois_v)
})
check_same <- try(print(pois_n), silent = TRUE)
test_that("Printing gives no error for type = none", {
  testthat::expect_identical(check_same, pois_n)
})
