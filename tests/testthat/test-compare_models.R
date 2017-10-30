context("compare_models")

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
larger <- adjust_loglik(gev_loglik, data = owtemps, init = init,
                        par_names = par_names)

# Restricted model
smaller_1 <- adjust_loglik(gev_loglik, data = owtemps, init = init,
                           par_names = par_names, fixed_pars = 6)

# Restricted model, using larger to start and character name
smaller_2 <- adjust_loglik(larger = larger, fixed_pars = "xi1")

res1 <- compare_models(larger, smaller_1)
res2 <- compare_models(larger, smaller_2)
res3 <- compare_models(larger, fixed_pars = 6)

my_tol <- 1e-5

test_that("approaches 1 and 3 agree", {
  testthat::expect_equal(res1, res3, tolerance = my_tol)
})
test_that("approaches 2 and 3 agree", {
  testthat::expect_equal(res2, res3, tolerance = my_tol)
})
