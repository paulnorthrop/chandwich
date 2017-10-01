context("compare_to_sandwich")

# Misspecified Poisson model for negative binomial data ----------
# ... following Section 5.1 of the
# "Object-Oriented Computation of Sandwich Estimators" vignette of the
# sandwich package
# https://cran.r-project.org/web/packages/sandwich/vignettes/sandwich-OOP.pdf

# Estimates and standard errors to reproduce
ests <- c(1.0633, 0.9961, -0.0491)
ses <- c(0.0414, 0.0535, 0.0231)
adj_ses <- c(0.0838, 0.1052, 0.0363)

set.seed(123)
x <- rnorm(250)
y <- rnbinom(250, mu = exp(1 + x), size = 1)
fm_pois <- glm(y ~ x + I(x^2), family = poisson)
adj_fn <- adjust_object(fm_pois)

pois_glm_loglik <- function(pars, y, x) {
  log_mu <- pars[1] + pars[2] * x + pars[3] * x ^ 2
  return(dpois(y, lambda = exp(log_mu), log = TRUE))
}
pois_res <- adjust_loglik(pois_glm_loglik, y = y, x = x, p = 3)

my_ests <- round(attr(pois_res, "MLE"), 4)
my_ses <- round(attr(pois_res, "SE"), 4)
my_adj_ses <- round(attr(pois_res, "adjSE"), 4)

my_tol <- 1e-5

test_that("MLEs agree", {
  testthat::expect_equal(ests, my_ests, tolerance = my_tol)
})
test_that("SEs agree", {
  testthat::expect_equal(ses, my_ses, tolerance = my_tol)
})
test_that("adjusted SEs agree", {
  testthat::expect_equal(adj_ses, my_adj_ses, tolerance = my_tol)
})
