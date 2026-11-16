## Threads inside a fit. The independent trials of a fit (starting values, candidate locations,
## spike branches) and the shapes of a profile run on ncores threads where the package is built
## with OpenMP, and must give bit-identical results for every number of threads; where it is not,
## ncores is accepted and has no effect, so the same comparisons hold. CRAN-mode tests use at most
## two threads.

same_fit <- function(a, b) {
  expect_identical(a$par, b$par)
  expect_identical(a$logLik, b$logLik)
  expect_identical(a$log_rates, b$log_rates)
  expect_identical(a$nfev, b$nfev)
  expect_identical(a$at_floor, b$at_floor)
  expect_identical(a$at_upper, b$at_upper)
  expect_identical(a$at_rate_bound, b$at_rate_bound)
  expect_identical(a$at_one_sided, b$at_one_sided)
}

test_that("a fit is identical on one and two threads", {
  set.seed(31)
  th <- rwsgn(12, mu = 1, beta = 1.6, lam = 0.3, sigma = 0.4)
  f1 <- fit_wsgn(th, se = "none", ncores = 1)
  f2 <- fit_wsgn(th, se = "none", ncores = 2)
  same_fit(f1, f2)
  expect_identical(f1$nfev, fit_wsgn(th, se = "none")$nfev)
  ## a profile too
  p1 <- profile_wsgn(th, c(0.8, 1.6), ncores = 1)
  p2 <- profile_wsgn(th, c(0.8, 1.6), ncores = 2)
  expect_identical(p1, p2)
})

test_that("ncores is validated, in fits and profiles", {
  th <- c(0.1, 0.5, 0.9, 1.2, 0.3, 0.7)
  expect_error(fit_wsgn(th, ncores = 0), "ncores")
  expect_error(fit_wsgn(th, ncores = 1.5), "ncores")
  expect_error(profile_wsgn(th, 1, ncores = 0), "ncores")
  expect_error(profile_wsgn(th, 1, ncores = NA), "ncores")
  expect_error(profile_wsgn(th, 1, ncores = c(1, 2)), "ncores")
})

test_that("fits are identical on one, two and four threads: exact, ties, grouped, one-sided", {
  skip_on_cran()
  set.seed(32)
  sets <- list(
    exact = list(th = rwsgn(60, mu = 1, beta = 1.5, lam = 0.35, sigma = 0.4)),
    low_shape = list(th = rwsgn(50, mu = 2, beta = 0.8, lam = 0.4, sigma = 0.3)),
    ties = list(th = round(rwsgn(70, mu = 1.5, beta = 1.2, lam = 0.5, sigma = 0.3), 2)),
    grouped = list(th = rwsgn(150, mu = 3, beta = 1.3, lam = 0.45, sigma = 0.5),
                   res = 2 * pi / 36),
    floor = list(th = rwsgn(40, mu = 1, beta = 0.6, lam = 0.3, sigma = 0.3), floor = 0.2)
  )
  x <- rwsgn(40, 1, 3, 0.5, sigma = 0.01)
  sets$one_sided <- list(th = signif(1 + abs(x - 1), 15))
  for (nm in names(sets)) {
    d <- sets[[nm]]
    fit <- function(nc) suppressMessages(
      fit_wsgn(d$th, beta_floor = d$floor, resolution = d$res, se = "none", ncores = nc))
    f1 <- fit(1)
    for (nc in c(2, 4)) {
      f <- fit(nc)
      same_fit(f1, f)
    }
  }
  th <- sets$low_shape$th
  p1 <- profile_wsgn(th, c(0.4, 0.8, 1.2, 2), ncores = 1)
  expect_identical(p1, profile_wsgn(th, c(0.4, 0.8, 1.2, 2), ncores = 2))
  expect_identical(p1, profile_wsgn(th, c(0.4, 0.8, 1.2, 2), ncores = 4))
})

test_that("the bootstrap does not depend on the number of threads", {
  skip_on_cran()
  set.seed(33)
  th <- rwsgn(20, mu = 1, beta = 1.6, lam = 0.3, sigma = 0.4)
  set.seed(1); b1 <- fit_wsgn(th, se = "bootstrap", B = 4, ncores = 1)
  set.seed(1); b2 <- fit_wsgn(th, se = "bootstrap", B = 4, ncores = 2)
  expect_identical(b1$boot, b2$boot)
  expect_identical(b1$se, b2$se)
})

## A likelihood evaluation with at least 128 observations (or 128 cells, when grouped) is computed in
## chunks of 64 as tasks and summed in index order on the calling thread; the value, and so the
## whole fit, must be bit-identical to the serial one.
test_that("chunked likelihood evaluations (n or cells at least 128) are identical on 1 and 2 threads", {
  set.seed(34)
  th <- rwsgn(130, mu = 1, beta = 1.6, lam = 0.3, sigma = 0.4)
  p1 <- profile_wsgn(th, 1.6, ncores = 1)
  expect_identical(p1, profile_wsgn(th, 1.6, ncores = 2))
  ## grouped: 180 cells of 2 degrees
  h <- 2 * pi / 180
  g <- (round(rwsgn(60, mu = 3, beta = 1.3, lam = 0.45, sigma = 0.5) / h) * h) %% (2 * pi)
  g1 <- suppressMessages(fit_wsgn(g, resolution = h, se = "none", ncores = 1))
  g2 <- suppressMessages(fit_wsgn(g, resolution = h, se = "none", ncores = 2))
  same_fit(g1, g2)
  pg <- profile_wsgn(g, 1.3, resolution = h, ncores = 1)
  expect_identical(pg, profile_wsgn(g, 1.3, resolution = h, ncores = 2))
})

test_that("chunked likelihood evaluations give identical fits on one, two and four threads", {
  skip_on_cran()
  set.seed(35)
  th <- rwsgn(140, mu = 1, beta = 1.6, lam = 0.3, sigma = 0.4)
  f1 <- suppressMessages(fit_wsgn(th, se = "none", ncores = 1))
  for (nc in c(2, 4)) same_fit(f1, suppressMessages(fit_wsgn(th, se = "none", ncores = nc)))
  h <- 2 * pi / 180
  g <- (round(rwsgn(100, mu = 3, beta = 1.3, lam = 0.45, sigma = 0.5) / h) * h) %% (2 * pi)
  g1 <- suppressMessages(fit_wsgn(g, resolution = h, se = "none", ncores = 1))
  for (nc in c(2, 4))
    same_fit(g1, suppressMessages(fit_wsgn(g, resolution = h, se = "none", ncores = nc)))
  p1 <- profile_wsgn(th, c(0.8, 1.6), ncores = 1)
  expect_identical(p1, profile_wsgn(th, c(0.8, 1.6), ncores = 4))
})

test_that("ncores above 256 is reduced with a message", {
  th <- c(0.1, 0.5, 0.9, 1.2, 0.3, 0.7)
  expect_message(p <- profile_wsgn(th, 1.5, ncores = 1e5), "256")
  expect_identical(p, profile_wsgn(th, 1.5))
})
