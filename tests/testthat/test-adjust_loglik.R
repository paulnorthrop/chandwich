context("adjust_loglik")

# ------------------------- Binomial model, rats data ----------------------

# Contributions to the independence loglikelihood
binom_loglik <- function(prob, data) {
  if (prob < 0 || prob > 1) {
    return(-Inf)
  }
  return(stats::dbinom(data[, "y"], data[, "n"], prob, log = TRUE))
}

# Different ways to set the number of parameters p

# Using par_names
rat_res_1 <- adjust_loglik(loglik = binom_loglik, data = rats, par_names = "p")

# Using init
rat_res_2 <- adjust_loglik(loglik = binom_loglik, data = rats, init = 0.1)

# Using p explicitly
rat_res_3 <- adjust_loglik(loglik = binom_loglik, data = rats, p = 1)

my_tol <- 1e-5

mle_1 <- as.numeric(attr(rat_res_1, "MLE"))
mle_2 <- attr(rat_res_2, "MLE")
mle_3 <- attr(rat_res_3, "MLE")
se_1 <- as.numeric(attr(rat_res_1, "SE"))
se_2 <- attr(rat_res_2, "SE")
se_3 <- attr(rat_res_3, "SE")
adjse_1 <- as.numeric(attr(rat_res_1, "adjSE"))
adjse_2 <- attr(rat_res_2, "adjSE")
adjse_3 <- attr(rat_res_3, "adjSE")

test_that("MLEs: par_names and init agree", {
  testthat::expect_equal(mle_1, mle_2, tolerance = my_tol)
})
test_that("MLEs: par_names and p agree", {
  testthat::expect_equal(mle_1, mle_3, tolerance = my_tol)
})
test_that("SEs: par_names and init agree", {
  testthat::expect_equal(se_1, se_2, tolerance = my_tol)
})
test_that("SEs: par_names and p agree", {
  testthat::expect_equal(se_1, se_3, tolerance = my_tol)
})
test_that("Adj. SEs: par_names and init agree", {
  testthat::expect_equal(adjse_1, adjse_2, tolerance = my_tol)
})
test_that("Adj. SEs: par_names and p agree", {
  testthat::expect_equal(adjse_1, adjse_3, tolerance = my_tol)
})

