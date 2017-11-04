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
conf_none <- conf_intervals(rat_res, type = "none")

check_same <- try(print(conf_v), silent = TRUE)
test_that("Printing gives no error for type = vertical", {
  testthat::expect_identical(check_same, conf_v)
})
check_same <- try(print(conf_none), silent = TRUE)
test_that("Printing gives no error for type = none", {
  testthat::expect_identical(check_same, conf_none)
})

