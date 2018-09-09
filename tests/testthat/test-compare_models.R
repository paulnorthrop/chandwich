context("compare_models and anova S3 method")

# GEV independence loglikelihood for the Oxford-Worthing annual maximum
# temperature dataset owtemps

gev_loglik <- function(pars, data) {
  o_pars <- pars[c(1, 3, 5)] + pars[c(2, 4, 6)]
  w_pars <- pars[c(1, 3, 5)] - pars[c(2, 4, 6)]
  if (o_pars[2] <= 0 | w_pars[2] <= 0) return(-Inf)
  o_data <- data[, "Oxford"]
  w_data <- data[, "Worthing"]
  check <- 1 + o_pars[3] * (o_data - o_pars[1]) / o_pars[2]
  if (any(check <= 0)) return(-Inf)
  check <- 1 + w_pars[3] * (w_data - w_pars[1]) / w_pars[2]
  if (any(check <= 0)) return(-Inf)
  o_loglik <- log_gev(o_data, o_pars[1], o_pars[2], o_pars[3])
  w_loglik <- log_gev(w_data, w_pars[1], w_pars[2], w_pars[3])
  return(o_loglik + w_loglik)
}

# Initial estimates (method of moments for the Gumbel case)
sigma <- as.numeric(sqrt(6 * diag(stats::var(owtemps))) / pi)
mu <- as.numeric(colMeans(owtemps) - 0.57722 * sigma)
init <- c(mean(mu), -diff(mu) / 2, mean(sigma), -diff(sigma) / 2, 0, 0)
par_names <-  c("mu0", "mu1", "sigma0", "sigma1", "xi0", "xi1")

# Full model
large <- adjust_loglik(gev_loglik, data = owtemps, init = init,
                        par_names = par_names)

# Restricted model, using larger to start and character name
medium <- adjust_loglik(larger = large, fixed_pars = "xi1")

# Restricted model
medium_2 <- adjust_loglik(gev_loglik, data = owtemps, init = init,
                          par_names = par_names, fixed_pars = 6)

# Different ways to make the same comparison
res1 <- compare_models(large, medium_2)
res2 <- compare_models(large, medium)
# Save for use later in anova.chandwich() tests
save_large_medium <- res2
res3 <- compare_models(large, fixed_pars = 6)
res4 <- compare_models(large, fixed_pars = "xi1")

my_tol <- 1e-5

test_that("approaches 1 and 3 agree", {
  testthat::expect_equal(res1, res3, tolerance = my_tol)
})
test_that("approaches 2 and 3 agree", {
  testthat::expect_equal(res2$p_value, res3$p_value, tolerance = my_tol)
})
test_that("approaches 3 and 4 agree", {
  testthat::expect_equal(res3$p_value, res4$p_value, tolerance = my_tol)
})
test_that("approaches 2 and 4 agree", {
  testthat::expect_equal(res2, res4, tolerance = my_tol)
})

# Repeat for approx = TRUE [approx only relevant when smaller is supplied]

medium_3 <- adjust_loglik(larger = large, fixed_pars = 6)

res1 <- compare_models(large, medium, approx = TRUE, method = "L-BFGS-B")
# Save for use later in anova.chandwich() tests
save_large_medium_approx <- res1
res2 <- compare_models(large, medium_3, approx = TRUE, method = "L-BFGS-B")

test_that("approx = TRUE, numeric and character which_pars agree", {
  testthat::expect_equal(res1, res2, tolerance = my_tol)
})

# Repeat fixed_pars numeric vs fixed_pars character, specifying an optim method

res3 <- compare_models(large, fixed_pars = 6, method = "L-BFGS-B")
res4 <- compare_models(large, fixed_pars = "xi1", method = "L-BFGS-B")

test_that("fixed_pars numeric vs character for L-BFGS-B", {
  testthat::expect_equal(res3$p_value, res4$p_value, tolerance = my_tol)
})

# Check that if type = "none" then approx = TRUE gives the same answer as
# approx = FALSE.  Also ask for "Nelder-Mead" (inappropriately when df = 1)
# and check that the code corrects for this

res1 <- compare_models(large, medium, approx = FALSE, type = "none",
                       method = "Nelder-Mead")
res2 <- compare_models(large, medium_3, approx = TRUE, type = "none",
                       method = "Nelder-Mead")

test_that("If type = none then approx has no effect", {
  testthat::expect_equal(res1, res2, tolerance = my_tol)
})

# Compare models where the larger model is not the full model

small <- adjust_loglik(larger = medium, fixed_pars = c("sigma1", "xi1"))
small_2 <- adjust_loglik(larger = medium, fixed_pars = c(4, 6))

res1 <- compare_models(medium, small)
# Save for use later in anova.chandwich() tests
save_medium_small <- res1
res2 <- compare_models(medium, small_2)

test_that("Larger model not full, numeric vs. character fixed_pars", {
  testthat::expect_equal(res1, res2, tolerance = my_tol)
})

# Repeat fixed_pars numeric vs fixed_pars character, specifying an optim method

small <- adjust_loglik(larger = medium, fixed_pars = c("sigma1", "xi1"),
                       method = "L-BFGS-B")
small_2 <- adjust_loglik(larger = medium, fixed_pars = c(4, 6),
                         method = "L-BFGS-B")

res1 <- compare_models(medium, small)
res2 <- compare_models(medium, small_2)

test_that("Larger model not full, numeric vs. character fixed_pars, L-BFGS-B", {
  testthat::expect_equal(res1, res2, tolerance = my_tol)
})

# Check that errors are produced when they should be

# Models are provided in the wrong order

check_error <- try(compare_models(medium, large), silent = TRUE)
test_that("Models are provided in the wrong order", {
  testthat::expect_identical(class(check_error), "try-error")
})

# Two models of the same size

medium_sigma1 <- adjust_loglik(larger = large, fixed_pars = "sigma1")
check_error <- try(compare_models(medium, medium_sigma1), silent = TRUE)
test_that("Models to be compared have the same size", {
  testthat::expect_identical(class(check_error), "try-error")
})

# Parameters fixed at different values

small_01 <- adjust_loglik(larger = medium, fixed_pars = c("sigma1", "xi1"),
                          fixed_at = c(0, 0.1))
check_error <- try(compare_models(medium, small_01), silent = TRUE)
test_that("Parameters fixed at different values", {
  testthat::expect_identical(class(check_error), "try-error")
})

# Values fixed in larger but not in smaller

medium_mu1 <- adjust_loglik(larger = large, fixed_pars = "mu1")
check_error <- try(compare_models(medium_mu1, small), silent = TRUE)
test_that("Values fixed in larger but not in smaller", {
  testthat::expect_identical(class(check_error), "try-error")
})

# Print

check_same <- utils::capture.output(print(res1))
check_res1 <- utils::capture.output(res1)
test_that("gev: print OK for character fixed_pars", {
  testthat::expect_identical(check_same, check_res1)
})

check_same <- utils::capture.output(print(res2))
check_res2 <- utils::capture.output(res2)
test_that("gev: print OK for numeric fixed_pars", {
  testthat::expect_identical(check_same, check_res2)
})

# Check that compare_models() and anova.chandwich() give the same results

tiny <- adjust_loglik(larger = small,
                      fixed_pars = c("mu1", "sigma1", "xi1"))
small_tiny <- compare_models(small, tiny)

# approx = FALSE

anova_res <- anova(large, medium, small, tiny)

test_that("anova ALRTS: large vs. medium", {
  testthat::expect_equal(save_large_medium$alrts, anova_res$ALRTS[2])
})
test_that("anova p-value: large vs. medium", {
  testthat::expect_equal(save_large_medium$p_value, anova_res$"Pr(>ALRTS)"[2])
})
test_that("anova ALRTS: medium vs. small", {
  testthat::expect_equal(save_medium_small$alrts, anova_res$ALRTS[3])
})
test_that("anova p-value: medium vs. small", {
  testthat::expect_equal(save_medium_small$p_value, anova_res$"Pr(>ALRTS)"[3])
})
test_that("anova ALRTS: small vs. tiny", {
  testthat::expect_equal(small_tiny$alrts, anova_res$ALRTS[4])
})
test_that("anova p-value: small vs. tiny", {
  testthat::expect_equal(small_tiny$p_value, anova_res$"Pr(>ALRTS)"[4])
})

# approx = TRUE

anova_res_approx <- anova(large, medium, small, tiny, approx = TRUE)
save_medium_small_approx <- compare_models(medium, small, approx = TRUE,
                                           method = "L-BFGS-B")
small_tiny_approx <- compare_models(small, tiny, approx = TRUE)

test_that("anova ALRTS: large vs. medium, approx", {
  testthat::expect_equal(save_large_medium_approx$alrts,
                         anova_res_approx$ALRTS[2])
})
test_that("anova p-value: large vs. medium, approx", {
  testthat::expect_equal(save_large_medium_approx$p_value,
                         anova_res_approx$"Pr(>ALRTS)"[2])
})
test_that("anova ALRTS: medium vs. small, approx", {
  testthat::expect_equal(save_medium_small_approx$alrts,
                         anova_res_approx$ALRTS[3])
})
test_that("anova p-value: medium vs. small, approx", {
  testthat::expect_equal(save_medium_small_approx$p_value,
                         anova_res_approx$"Pr(>ALRTS)"[3])
})
test_that("anova ALRTS: medium vs. small, approx", {
  testthat::expect_equal(small_tiny_approx$alrts,
                         anova_res_approx$ALRTS[4])
})
test_that("anova p-value: medium vs. small, approx", {
  testthat::expect_equal(small_tiny_approx$p_value,
                         anova_res_approx$"Pr(>ALRTS)"[4])
})

