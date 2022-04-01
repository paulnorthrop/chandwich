#context("plot.chandwich")

# ------------------------- Binomial model, rats data ----------------------

# Contributions to the independence loglikelihood
binom_loglik <- function(prob, data) {
  if (prob < 0 || prob > 1) {
    return(-Inf)
  }
  return(dbinom(data[, "y"], data[, "n"], prob, log = TRUE))
}
rat_res <- adjust_loglik(loglik = binom_loglik, data = rats, par_names = "p")

# Vertically adjusted loglikelihood only
check_NULL <- try(plot(rat_res), silent = TRUE)
test_that("Default plot gives no error", {
  testthat::expect_identical(check_NULL, NULL)
})
check_NULL <- try(plot(rat_res, type = 1:4, legend_pos = "bottom"),
                  silent = TRUE)
test_that("Adding arguments gives no error", {
  testthat::expect_identical(check_NULL, NULL)
})
check_NULL <- try(plot(rat_res, type = 1:4, xlim = c(0, 1),
                       legend_pos = "bottom"), silent = TRUE)
test_that("Adding arguments, inc. xlim, gives no error", {
  testthat::expect_identical(check_NULL, NULL)
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

# Quadratic model
pois_quad <- adjust_loglik(pois_glm_loglik, y = y, x = x, p = 3)
check_error <- try(plot(pois_quad), silent = TRUE)

test_that("More than one free parameter produces an error", {
  testthat::expect_identical(class(check_error), "try-error")
})
