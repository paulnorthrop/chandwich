context("compare_to_sandwich")

# Misspecified Poisson model for negative binomial data ----------
# ... following Section 5.1 of the
# "Object-Oriented Computation of Sandwich Estimators" vignette of the
# sandwich package
# https://cran.r-project.org/web/packages/sandwich/vignettes/sandwich-OOP.pdf

my_tol <- 1e-5

# Estimates and standard errors to reproduce
ests <- c(1.0633, 0.9961, -0.0491)
ses <- c(0.0414, 0.0535, 0.0231)
adj_ses <- c(0.0838, 0.1052, 0.0363)

set.seed(123)
x <- stats::rnorm(250)
y <- stats::rnbinom(250, mu = exp(1 + x), size = 1)
fm_pois <- stats::glm(y ~ x + I(x^2), family = poisson)

pois_glm_loglik <- function(pars, y, x) {
  log_mu <- pars[1] + pars[2] * x + pars[3] * x ^ 2
  return(stats::dpois(y, lambda = exp(log_mu), log = TRUE))
}
pois_res <- adjust_loglik(pois_glm_loglik, y = y, x = x, p = 3)
my_ests <- round(attr(pois_res, "MLE"), 4)
my_ses <- round(attr(pois_res, "SE"), 4)
my_adj_ses <- round(attr(pois_res, "adjSE"), 4)

test_that("MLEs agree", {
  testthat::expect_equal(ests, my_ests, tolerance = my_tol)
})
test_that("SEs agree", {
  testthat::expect_equal(ses, my_ses, tolerance = my_tol)
})
test_that("adjusted SEs agree", {
  testthat::expect_equal(adj_ses, my_adj_ses, tolerance = my_tol)
})

# Repeat when supplying the MLE

pois_res_2 <- adjust_loglik(pois_glm_loglik, y = y, x = x, p = 3,
                            mle = fm_pois$coefficients)
my_ests_2a <- round(attr(pois_res_2, "MLE"), 4)
my_ses_2a <- round(attr(pois_res_2, "SE"), 4)
my_adj_ses_2a <- round(attr(pois_res_2, "adjSE"), 4)

test_that("MLEs agree, mle and H supplied", {
  testthat::expect_equal(my_ests, my_ests_2a, tolerance = my_tol)
})
test_that("SEs agree, mle and H supplied", {
  testthat::expect_equal(my_ses, my_ses_2a, tolerance = my_tol)
})
test_that("adjusted SEs agree, mle and H supplied", {
  testthat::expect_equal(my_adj_ses, my_adj_ses_2a, tolerance = my_tol)
})

# Repeat when supplying the MLE and the Hessian H of the loglikelihood
# evaluated at the MLE

pois_res_2 <- adjust_loglik(pois_glm_loglik, y = y, x = x, p = 3,
                            mle = fm_pois$coefficients,
                            H = -solve(stats::vcov(fm_pois)))
my_ests_2 <- round(attr(pois_res_2, "MLE"), 4)
my_ses_2 <- round(attr(pois_res_2, "SE"), 4)
my_adj_ses_2 <- round(attr(pois_res_2, "adjSE"), 4)

test_that("MLEs agree, mle and H supplied", {
  testthat::expect_equal(my_ests, my_ests_2, tolerance = my_tol)
})
test_that("SEs agree, mle and H supplied", {
  testthat::expect_equal(my_ses, my_ses_2, tolerance = my_tol)
})
test_that("adjusted SEs agree, mle and H supplied", {
  testthat::expect_equal(my_adj_ses, my_adj_ses_2, tolerance = my_tol)
})

# Repeat when supplying the MLE and the Hessian H of the loglikelihood
# evaluated at the MLE, and the variance V of the vector of cluster-specific
# contributions to the score vector.
#
# Use sandwich::bread() and sandwich::meat() to evaluate H and V respectively.

got_sandwich <- requireNamespace("sandwich", quietly = TRUE)
print(got_sandwich)
if (got_sandwich) {
  n_obs <- stats::nobs(fm_pois)
  pois_res_2b <- adjust_loglik(pois_glm_loglik, y = y, x = x, p = 3,
                               mle = fm_pois$coefficients,
                               H = -solve(sandwich::bread(fm_pois) / n_obs),
                               V = sandwich::meat(fm_pois) * n_obs)
  my_ests_2b <- round(attr(pois_res_2b, "MLE"), 4)
  my_ses_2b <- round(attr(pois_res_2b, "SE"), 4)
  my_adj_ses_2b <- round(attr(pois_res_2b, "adjSE"), 4)

  test_that("MLEs agree, mle and H supplied", {
    testthat::expect_equal(my_ests, my_ests_2b, tolerance = my_tol)
  })
  test_that("SEs agree, mle and H supplied", {
    testthat::expect_equal(my_ses, my_ses_2b, tolerance = my_tol)
  })
  test_that("adjusted SEs agree, mle and H supplied", {
    testthat::expect_equal(my_adj_ses, my_adj_ses_2b, tolerance = my_tol)
  })
}

# Repeat using algebraic derivatives and Hessian

pois_alg_deriv <- function(pars, y, x) {
  mu <- exp(pars[1] + pars[2] * x + pars[3] * x ^ 2)
  return(cbind(y - mu, x * (y - mu), x ^ 2 * (y - mu)))
}

pois_alg_hess <- function(pars, y, x) {
  mu <- exp(pars[1] + pars[2] * x + pars[3] * x ^ 2)
  alg_hess <- matrix(0, 3, 3)
  alg_hess[1, ] <- -c(sum(mu), sum(x * mu), sum(x ^ 2 * mu))
  alg_hess[2, ] <- -c(sum(x * mu), sum(x ^ 2 * mu), sum(x ^ 3 * mu))
  alg_hess[3, ] <- -c(sum(x ^ 2 * mu), sum(x ^ 3 * mu), sum(x ^ 4 * mu))
  return(alg_hess)
}

pois_res_alg <- adjust_loglik(pois_glm_loglik, y = y, x = x, p = 3,
                              alg_deriv = pois_alg_deriv,
                              alg_hess = pois_alg_hess)

my_ests_alg <- round(attr(pois_res_alg, "MLE"), 4)
my_ses_alg <- round(attr(pois_res_alg, "SE"), 4)
my_adj_ses_alg <- round(attr(pois_res_alg, "adjSE"), 4)

test_that("MLEs agree, algebraic derivatives", {
  testthat::expect_equal(ests, my_ests_alg, tolerance = my_tol)
})
test_that("SEs agree, algebraic derivatives", {
  testthat::expect_equal(ses, my_ses_alg, tolerance = my_tol)
})
test_that("adjusted SEs agree, algebraic derivatives", {
  testthat::expect_equal(adj_ses, my_adj_ses_alg, tolerance = my_tol)
})

# Check that for the linear model we get the same answer using algebraic
# derivatives as using numerical derivatives

pois_lin <- adjust_loglik(larger = pois_res, fixed_pars = 3)
pois_lin_alg <- adjust_loglik(larger = pois_res_alg, fixed_pars = 3)

lin_ests <- round(attr(pois_lin, "MLE"), 4)
lin_ses <- round(attr(pois_lin, "SE"), 4)
lin_adj_ses <- round(attr(pois_lin, "adjSE"), 4)
lin_ests_alg <- round(attr(pois_lin, "MLE"), 4)
lin_ses_alg <- round(attr(pois_lin, "SE"), 4)
lin_adj_ses_alg <- round(attr(pois_lin, "adjSE"), 4)

test_that("MLEs agree: linear, algebraic vs numeric derivatives", {
  testthat::expect_equal(lin_ests, lin_ests_alg, tolerance = my_tol)
})
test_that("SEs agree: linear, algebraic vs numeric derivatives", {
  testthat::expect_equal(lin_ses, lin_ses_alg, tolerance = my_tol)
})
test_that("adj. SEs agree: linear, algebraic vs numeric derivatives", {
  testthat::expect_equal(lin_adj_ses, lin_adj_ses_alg, tolerance = my_tol)
})
