context("conf_reg")

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
pois_vertical <- conf_region(pois_lin)
pois_none <- conf_region(pois_lin, type = "none")
check_NULL <- try(plot(pois_none, pois_vertical, conf = c(50, 75, 95, 99),
                       col = 2:1, lwd = 2, lty = 1), silent = TRUE)


test_that("Consistent names gives no error", {
  testthat::expect_identical(check_NULL, NULL)
})

new_pois_lin <- adjust_loglik(pois_glm_loglik, y = y, x = x, par_names = pars,
                              fixed_pars = "gamma", name = "wrong_name")
new_pois_none <- conf_region(new_pois_lin, type = "none")
check_error <- try(plot(new_pois_none, pois_vertical, conf = c(50, 75, 95, 99)),
                   silent = TRUE)

test_that("Inconsistent names gives an error", {
  testthat::expect_identical(class(check_error), "try-error")
})
