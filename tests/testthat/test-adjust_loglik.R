context("adjust_loglik")

# ------------------------- Binomial model, rats data ----------------------

# Contributions to the independence loglikelihood
binom_loglik <- function(prob, data) {
  if (prob < 0 || prob > 1) {
    return(-Inf)
  }
  return(stats::dbinom(data[, "y"], data[, "n"], prob, log = TRUE))
}

# Different ways to set the number of parameters p and the cluster indicator

# Using par_names
rat_res_1 <- adjust_loglik(loglik = binom_loglik, data = rats, par_names = "p")

# Using init
rat_res_2 <- adjust_loglik(loglik = binom_loglik, data = rats, init = 0.1)

# Using p explicitly
rat_res_3 <- adjust_loglik(loglik = binom_loglik, data = rats, p = 1)

# Setting cluster explicitly
rat_res_4 <- adjust_loglik(loglik = binom_loglik, data = rats, p = 1,
                           cluster = 1:dim(rats)[1])

my_tol <- 1e-5

mle_1 <- as.numeric(attr(rat_res_1, "MLE"))
mle_2 <- attr(rat_res_2, "MLE")
mle_3 <- attr(rat_res_3, "MLE")
mle_4 <- attr(rat_res_4, "MLE")
se_1 <- as.numeric(attr(rat_res_1, "SE"))
se_2 <- attr(rat_res_2, "SE")
se_3 <- attr(rat_res_3, "SE")
se_4 <- attr(rat_res_4, "SE")
adjse_1 <- as.numeric(attr(rat_res_1, "adjSE"))
adjse_2 <- attr(rat_res_2, "adjSE")
adjse_3 <- attr(rat_res_3, "adjSE")
adjse_4 <- attr(rat_res_4, "adjSE")

test_that("MLEs: par_names and init agree", {
  testthat::expect_equal(mle_1, mle_2, tolerance = my_tol)
})
test_that("MLEs: par_names and p agree", {
  testthat::expect_equal(mle_1, mle_3, tolerance = my_tol)
})
test_that("MLEs: cluster and default agree", {
  testthat::expect_equal(mle_3, mle_4, tolerance = my_tol)
})
test_that("SEs: par_names and init agree", {
  testthat::expect_equal(se_1, se_2, tolerance = my_tol)
})
test_that("SEs: par_names and p agree", {
  testthat::expect_equal(se_1, se_3, tolerance = my_tol)
})
test_that("SEs: cluster and default agree", {
  testthat::expect_equal(se_3, se_4, tolerance = my_tol)
})
test_that("Adj. SEs: par_names and init agree", {
  testthat::expect_equal(adjse_1, adjse_2, tolerance = my_tol)
})
test_that("Adj. SEs: par_names and p agree", {
  testthat::expect_equal(adjse_1, adjse_3, tolerance = my_tol)
})
test_that("Adj. SEs: cluster and default agree", {
  testthat::expect_equal(adjse_3, adjse_4, tolerance = my_tol)
})

# Repeat 1 and 2 when using method = "L-BFGS-B"

# Different ways to set the number of parameters p and the cluster indicator
# Using par_names
rat_res_1 <- adjust_loglik(loglik = binom_loglik, data = rats, par_names = "p",
                           method = "L-BFGS-B")
# Using init
rat_res_2 <- adjust_loglik(loglik = binom_loglik, data = rats, init = 0.1,
                           method = "L-BFGS-B")

mle_1 <- as.numeric(attr(rat_res_1, "MLE"))
mle_2 <- attr(rat_res_2, "MLE")
se_1 <- as.numeric(attr(rat_res_1, "SE"))
se_2 <- attr(rat_res_2, "SE")
adjse_1 <- as.numeric(attr(rat_res_1, "adjSE"))
adjse_2 <- attr(rat_res_2, "adjSE")
test_that("MLEs: par_names and init agree, method = L-BFGS-B", {
  testthat::expect_equal(mle_1, mle_2, tolerance = my_tol)
})
test_that("SEs: par_names and init agree, method = L-BFGS-B", {
  testthat::expect_equal(se_1, se_2, tolerance = my_tol)
})
test_that("Adj. SEs: par_names and init agree, method = L-BFGS-B", {
  testthat::expect_equal(adjse_1, adjse_2, tolerance = my_tol)
})

# Check that errors are produced when they should be

# Lengths of par_names and init are not consistent
check_error <- try(adjust_loglik(loglik = binom_loglik, data = rats,
                                 par_names = "p", init = c(0.1, 0.1)),
                                 silent = TRUE)
test_that("Lengths of par_names and init are not consistent", {
  testthat::expect_identical(class(check_error), "try-error")
})

# p and length of par_names are not consistent
check_error <- try(adjust_loglik(loglik = binom_loglik, data = rats, p = 2,
                                 par_names = "p"),
                   silent = TRUE)
test_that("p and length of par_names are not consistent", {
  testthat::expect_identical(class(check_error), "try-error")
})

# p and length of init are not consistent
check_error <- try(adjust_loglik(loglik = binom_loglik, data = rats, p = 1,
                                 init = c(0.1, 0.1)),
                   silent = TRUE)
test_that("p and length of init are not consistent", {
  testthat::expect_identical(class(check_error), "try-error")
})

# p and lengths of par_names and init are not consistent
check_error <- try(adjust_loglik(loglik = binom_loglik, data = rats, p = 2,
                                 par_names = "p", init = 0.1),
                   silent = TRUE)
test_that("p, lengths of par_names and init are not consistent", {
  testthat::expect_identical(class(check_error), "try-error")
})
