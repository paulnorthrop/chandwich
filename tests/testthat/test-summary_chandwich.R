context("summary_chandwich")

# ------------------------- Binomial model, rats data ----------------------

# Contributions to the independence loglikelihood
binom_loglik <- function(prob, data) {
  if (prob < 0 || prob > 1) {
    return(-Inf)
  }
  return(dbinom(data[, "y"], data[, "n"], prob, log = TRUE))
}

rat_res_1 <- adjust_loglik(loglik = binom_loglik, data = rats, par_names = "p")
res1 <- as.numeric(summary(rat_res_1))

rat_res_2 <- adjust_loglik(loglik = binom_loglik, data = rats, p = 1)
res2 <- as.numeric(summary(rat_res_2))

rat_res_3 <- adjust_loglik(loglik = binom_loglik, data = rats, init = 0.1)
res3 <- as.numeric(summary(rat_res_3))

test_that("par_names and p are equivalent", {
  testthat::expect_identical(res1, res2)
})
test_that("p and init are equivalent", {
  testthat::expect_identical(res2, res3)
})
