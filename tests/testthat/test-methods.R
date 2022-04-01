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
pars <- c("alpha", "beta", "gamma")
pois_quad <- adjust_loglik(pois_glm_loglik, y = y, x = x, p = 3,
                           par_names = pars)
# Linear model (gamma fixed at 0).
pois_lin <- adjust_loglik(pois_glm_loglik, y = y, x = x, par_names = pars,
                          fixed_pars = "gamma")

# ==================================== coef ===================================

#context("coef.chandwich")

test_that("coef.chandwich: complete = TRUE", {
  testthat::expect_identical(attributes(pois_lin)$res_MLE,
                         coef(pois_lin, complete = TRUE))
})
test_that("coef.chandwich: complete = FALSE", {
  testthat::expect_identical(attributes(pois_lin)$MLE,
                         coef(pois_lin))
})

# ==================================== vcov ===================================

#context("vcov.chandwich")

test_that("vcov.chandwich: adjusted = TRUE", {
  testthat::expect_identical(attributes(pois_lin)$adjVC,
                             vcov(pois_lin))
})
test_that("vcov.chandwich: adjusted = FALSE", {
  testthat::expect_identical(attributes(pois_lin)$VC,
                             vcov(pois_lin, adjusted = FALSE))
})
test_that("vcov.chandwich: complete = TRUE, column", {
  testthat::expect_equal(vcov(pois_lin, complete = TRUE)[, "gamma"],
                              as.numeric(rep(NA, 3)), ignore_attr = TRUE)
})
test_that("vcov.chandwich: complete = TRUE, row", {
  testthat::expect_equal(vcov(pois_lin, complete = TRUE)["gamma", ],
                              as.numeric(rep(NA, 3)), ignore_attr = TRUE)
})

# =================================== logLik ==================================

#context("logLik.chandwich")

test_that("logLik.chandwich: logLik and loglik at MLE agree", {
  testthat::expect_identical(as.numeric(logLik(pois_lin)),
                             pois_lin(attr(pois_lin, "MLE")))
})
test_that("logLik.chandwich: quadratic logLik > linear logLik", {
  testthat::expect_gte(logLik(pois_quad), logLik(pois_lin))
})
test_that("logLik.chandwich: quadratic df", {
  testthat::expect_equal(attr(logLik(pois_quad), "df"), 3)
})
test_that("logLik.chandwich: linear df", {
  testthat::expect_equal(attr(logLik(pois_lin), "df"), 2)
})


