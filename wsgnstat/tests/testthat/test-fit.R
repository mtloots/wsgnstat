## Maximum likelihood for the four-parameter wrapped distribution: the shape floor, the location
## search, the grouped likelihood, information and bootstrap standard errors. Where the circular
## package is available the fits are checked against values computed independently (in Python, by
## the authors' own code) for the periwinkle and ant data of Fisher (1993).

test_that("a fit with shape above one recovers the truth, with information standard errors", {
  set.seed(5)
  truth <- c(mu = 1, sigma = 0.5, lam = 0.3, beta = 1.8)
  th <- rwsgn(100, truth["mu"], truth["beta"], truth["lam"], sigma = truth["sigma"])
  f <- fit_wsgn(th)
  expect_s3_class(f, "wsgn_fit")
  expect_identical(f$se_method, "information")
  expect_named(coef(f), c("mu", "sigma", "lam", "beta"))
  expect_true(all(abs(coef(f) - truth) < 4 * f$se))
  expect_false(f$at_floor)
  cf <- coef(f)
  expect_equal(unname(f$rates), unname(c(cf["lam"], 1 - cf["lam"]) / cf["sigma"]),
               tolerance = 1e-14)
  ## the log-likelihood is the sum of the log-densities at the estimate, and it beats the truth
  ll <- sum(dwsgn(th, cf["mu"], cf["beta"], cf["lam"], sigma = cf["sigma"], log = TRUE))
  expect_equal(f$logLik, ll, tolerance = 1e-12)
  expect_gte(f$logLik, sum(dwsgn(th, 1, 1.8, 0.3, sigma = 0.5, log = TRUE)))
  expect_identical(as.numeric(logLik(f)), f$logLik)
  expect_identical(attr(logLik(f), "df"), 4L)
  expect_identical(dim(vcov(f)), c(4L, 4L))
  expect_equal(sqrt(diag(vcov(f))), f$se, tolerance = 1e-14, ignore_attr = TRUE)
  expect_output(print(f), "information")
})

test_that("the expected information agrees with an independent computation", {
  set.seed(5)
  f <- fit_wsgn(rwsgn(100, 1, 1.8, 0.3, sigma = 0.5))
  par <- coef(f)
  th <- par[1] - pi + (seq_len(4000) - 0.5) * 2 * pi / 4000
  lg <- function(p) dwsgn(th, p[1], p[4], p[3], sigma = p[2], log = TRUE)
  sc <- vapply(1:4, function(j) {
    h <- 1e-4
    e <- replace(numeric(4), j, h)
    (lg(par + e) - lg(par - e)) / (2 * h)
  }, numeric(length(th)))
  info <- crossprod(sc * sqrt(exp(lg(par)) * 2 * pi / 4000))
  ## vcov(f) is I(par)^-1 / n
  expect_lt(max(abs(solve(vcov(f)) / 100 - info) / abs(info)), 1e-3)
  expect_equal(unname(f$se), sqrt(diag(solve(info)) / 100), tolerance = 1e-3)
})

test_that("at a fitted shape below one the location is an observation and the bootstrap is used", {
  skip_on_cran()
  ## a deterministic sample: quantiles of van der Corput points (as in the Python tests)
  vdc <- function(n) rowSums(outer(seq_len(n), 1:20, function(i, j) ((i %/% 2^(j - 1)) %% 2) / 2^j))
  th <- qwsgn(vdc(90), 2, 0.6, 0.4, sigma = 0.25)
  set.seed(8)
  expect_message(fit_wsgn(th), "at most one")
  f <- fit_wsgn(th, se = "bootstrap", B = 4)
  expect_lte(coef(f)["beta"], 1)
  expect_true(coef(f)["mu"] %in% (th %% (2 * pi)))
  expect_identical(f$se_method, "bootstrap")
  expect_identical(dim(f$boot), c(4L, 5L))
  expect_true(all(is.finite(f$se) & f$se >= 0)); expect_true(all(f$se[c("mu", "sigma", "lam")] > 0))
  expect_gte(f$frac_beta_lt1, 0)
  expect_lte(f$frac_beta_lt1, 1)
  expect_null(vcov(f))
  ## the same bootstrap samples refitted with two threads give the same refits
  set.seed(9)
  b1 <- fit_wsgn(th, se = "bootstrap", B = 2, ncores = 1)$boot
  set.seed(9)
  b2 <- fit_wsgn(th, se = "bootstrap", B = 2, ncores = 2)$boot
  expect_identical(b1, b2)
})

test_that("a fit on the shape floor is flagged, and the symmetric member fixes lam", {
  set.seed(10)
  th <- rwsgn(80, mu = 1, beta = 0.5, lam = 0.5, sigma = 0.05)
  f <- fit_wsgn(th, beta_floor = 1.5, se = "none")
  expect_true(f$at_floor)
  expect_equal(unname(coef(f)["beta"]), 1.5, tolerance = 1e-5)
  expect_output(print(f), "ON the floor")
  fs <- fit_wsgn(rwsgn(60, 1, 2, 0.5, sigma = 0.5), symmetric = TRUE)
  expect_identical(unname(coef(fs)["lam"]), 0.5)
  expect_true(is.na(fs$se["lam"]))
  expect_identical(attr(logLik(fs), "df"), 3L)
})

test_that("the profile likelihood is the fit at fixed shapes and rises without bound below", {
  skip_on_cran()
  set.seed(11)
  th <- rwsgn(60, mu = 1, beta = 1.2, lam = 0.35, sigma = 0.4)
  f <- fit_wsgn(th, se = "none")
  pr <- profile_wsgn(th, beta = c(0.02, unname(coef(f)["beta"]), 3))
  expect_equal(pr$logLik[2], f$logLik, tolerance = 1e-6)
  expect_gt(pr$logLik[1], f$logLik)
  expect_lt(pr$logLik[3], f$logLik)
  expect_error(profile_wsgn(th, -1), "positive")
})

test_that("the grouped likelihood is the likelihood of the cells, and its fit is bounded", {
  set.seed(12)
  h <- pi / 18
  th <- (round(rwsgn(120, mu = 3, beta = 1.2, lam = 0.4, sigma = 0.3) / h) * h) %% (2 * pi)
  ll <- grouped_loglik(th, h, mu = 3, sigma = 0.3, lam = 0.4, beta = 1.2)
  cell <- round(th / h) %% 36
  pr <- pwsgn(((0:35) + 0.5) * h, 3, 1.2, 0.4, sigma = 0.3) -
    pwsgn(((0:35) - 0.5) * h, 3, 1.2, 0.4, sigma = 0.3)
  pr[1] <- pr[1] + 1   # the cell centred at 0 wraps: [-h/2, 0) is [2 pi - h/2, 2 pi)
  pr[1] <- pr[1] - pwsgn(2 * pi - h / 2, 3, 1.2, 0.4, sigma = 0.3)
  expect_equal(ll, sum(table(factor(cell, 0:35)) * log(pr)), tolerance = 1e-10)
  g <- fit_wsgn(th, resolution = h, se = "none")
  expect_true(g$grouped)
  expect_lte(g$logLik, 0)
  expect_gte(g$logLik, ll)
  expect_equal(g$logLik, grouped_loglik(th, h, coef(g)["mu"], coef(g)["sigma"], coef(g)["lam"],
                                        coef(g)["beta"]), tolerance = 1e-12)
  expect_error(fit_wsgn(th, resolution = 0.7), "whole number")
  expect_error(grouped_loglik(th, 4, 3, 0.3, 0.4, 1.2), "resolution")
})

test_that("grouped bootstrap standard errors refit rounded samples", {
  skip_on_cran()
  set.seed(13)
  h <- pi / 6
  th <- (round(rwsgn(50, mu = 1, beta = 1.5, lam = 0.5, sigma = 0.5) / h) * h) %% (2 * pi)
  g <- fit_wsgn(th, resolution = h, se = "bootstrap", B = 3)
  expect_identical(g$se_method, "bootstrap")
  expect_identical(nrow(g$boot), 3L)
  expect_true(all(is.finite(g$se)))
})

test_that("the periwinkle and ant fits agree with independent computations", {
  skip_on_cran()
  skip_if_not_installed("circular")
  e <- new.env()
  utils::data("fisherB20", "fisherB7", package = "circular", envir = e)
  peri <- e$fisherB20$theta * pi / 180
  ants <- as.numeric(e$fisherB7) * pi / 180
  ## exact-data fit under the floor 0.3: the location lands on the double tie at 60 degrees
  f <- fit_wsgn(peri, beta_floor = 0.3, se = "none")
  expect_identical(unname(coef(f)["mu"]), 60 * pi / 180)
  expect_equal(unname(coef(f)), c(1.0471975511965979, 0.0810581284965923, 0.12834861700500036,
                                  0.8730181786484532), tolerance = 1e-5)
  expect_equal(f$logLik, -27.296644722522444, tolerance = 1e-8)
  ## profile likelihood of the periwinkles: the spike branch overtakes below shape 0.15
  pp <- profile_wsgn(peri, c(1, 0.3, 0.05))
  expect_equal(pp$logLik, c(-27.343, -28.396, 4.185), tolerance = 1e-4)
  ## grouped fit of the ants (10-degree cells): an interior maximum
  g <- fit_wsgn(ants, resolution = pi / 18, se = "none")
  expect_equal(unname(coef(g)), c(3.1921560041866837, 0.08321258372285376, 0.48184120043855616,
                                  0.5085743634327373), tolerance = 1e-5)
  expect_equal(g$logLik, -306.9329052056407, tolerance = 1e-9)
  expect_false(g$at_floor)
})

test_that("the fit is the maximum over the floor branch as well", {
  skip_on_cran()
  skip_if_not_installed("circular")
  e <- new.env()
  utils::data("fisherB20", package = "circular", envir = e)
  peri <- e$fisherB20$theta * pi / 180
  ## at low floors the spike on a tie beats the regular fit: the fit is on the floor and equals the
  ## profile there (the triple tie at 98 degrees at 0.10, the double tie at 60 degrees at 0.15)
  for (fl in c(0.10, 0.15)) {
    f <- fit_wsgn(peri, beta_floor = fl, se = "none")
    expect_true(f$at_floor)
    expect_equal(f$logLik, profile_wsgn(peri, fl)$logLik, tolerance = 1e-8)
  }
  f10 <- fit_wsgn(peri, beta_floor = 0.10, se = "none")
  expect_equal(f10$logLik, -21.007222, tolerance = 1e-6)
  expect_identical(unname(coef(f10)["mu"]), 98 * pi / 180)
  expect_equal(fit_wsgn(peri, beta_floor = 0.15, se = "none")$logLik, -27.010481, tolerance = 1e-6)
  ## at 0.3 the regular fit wins and is unchanged
  f30 <- fit_wsgn(peri, beta_floor = 0.3, se = "none")
  expect_false(f30$at_floor)
  expect_equal(f30$logLik, -27.296644722522444, tolerance = 1e-8)
  expect_identical(unname(coef(f30)["mu"]), 60 * pi / 180)
  ## a fit on a bound gets bootstrap, never information, standard errors
  expect_error(fit_wsgn(peri, beta_floor = 0.10, se = "information"), "bound")
})

test_that("the search reaches the one-sided branch and the profile equals the floor search", {
  skip_on_cran()
  ## a one-sided maximum at shape about 0.83 that the regular starts used to miss
  set.seed(460)
  th <- rwsgn(60, 1, 1, 0.2, sigma = 0.3)
  f <- suppressMessages(fit_wsgn(th))
  expect_gte(f$logLik, profile_wsgn(th, 0.8)$logLik - 1e-8)
  expect_gte(f$logLik, -87.0644)
  ## a spike on a near-tie: the profile at the floor equals the fit's floor search
  set.seed(7017)
  th <- rwsgn(60, 1, 0.8, 0.3, sigma = 0.7)
  f <- suppressMessages(fit_wsgn(th))
  expect_true(f$at_floor)
  expect_equal(profile_wsgn(th, 0.1)$logLik, f$logLik, tolerance = 1e-8)
})

test_that("bootstrap draws come from the rates, so a degenerate fit can be bootstrapped", {
  skip_on_cran()
  set.seed(7021)
  th <- rwsgn(60, 1, 0.8, 0.5, sigma = 0.7)
  expect_message(f0 <- fit_wsgn(th), "bound")
  expect_null(f0$se)
  expect_true(all(is.finite(f0$log_rates)))
  expect_equal(unname(f0$rates), exp(unname(f0$log_rates)))
  f <- fit_wsgn(th, se = "bootstrap", B = 2)
  expect_identical(nrow(f$boot), 2L)
  expect_true(all(is.finite(f$boot)))
})

test_that("a concentrated fit keeps finite information standard errors", {
  set.seed(1)
  th <- rwsgn(50, 1, 2, 0.3, sigma = 5e-6)
  f <- fit_wsgn(th, ngrid = 2000)
  expect_identical(f$se_method, "information")
  expect_true(all(is.finite(f$se)))
  ## the standard error of the location scales with sigma
  expect_lt(f$se["mu"] / coef(f)["sigma"], 5)
})

test_that("invalid fitting arguments are refused", {
  expect_error(fit_wsgn(rep(2, 10)), "three distinct")
  expect_error(fit_wsgn(c(1, 1, 2, 2)), "three distinct")
  expect_error(fit_wsgn(c(0.1, 0.5, Inf, 2)), "finite")
  expect_error(fit_wsgn(c(0.1, 0.5, NaN, 2)), "finite")
  expect_error(fit_wsgn(rep(c(1, 1.01), 5), resolution = pi / 18), "single cell")
  expect_error(grouped_loglik(c(0.1, Inf, 1), pi / 18, 1, 0.4, 0.3, 0.5), "finite")
  expect_error(grouped_loglik(c(1, 2, 3), NULL, 1, 0.4, 0.3, 1.5), "resolution")
  expect_error(fit_wsgn(1:10, K = 1e10), "K")
  expect_error(fit_wsgn(1:10, beta_floor = 25), "beta_floor")
  expect_error(fit_wsgn(1:10, se = "bootstrap", B = NA), "B")
  expect_error(fit_wsgn(1:10, window = NA), "window")
  expect_error(fit_wsgn(1:10, ngrid = 0), "ngrid")
  expect_error(fit_wsgn(c(1, 2)), "three")
  expect_error(fit_wsgn(c(1, 2, NA)), "three")
  expect_error(fit_wsgn(1:10, beta_floor = -1), "beta_floor")
  expect_error(fit_wsgn(1:10, K = 0), "K")
  expect_error(fit_wsgn(1:10, se = "bootstrap", B = 1), "B")
})

test_that("a one-sided fit short of the rate bound is flagged and gets no information errors", {
  skip_on_cran()
  for (sd in c(3, 4, 20)) {
    set.seed(sd)
    th <- rwsgn(100, 1, 3, 0.4, sigma = 0.01)
    f <- suppressMessages(fit_wsgn(th))
    expect_true(f$at_one_sided)
    expect_false(f$at_rate_bound)
    expect_gt(abs(diff(f$log_rates)), log(1e12))
    expect_identical(f$se_method, "none")
    expect_null(f$se)
    expect_output(print(f), "one-sided")
  }
  set.seed(3)
  th <- rwsgn(100, 1, 3, 0.4, sigma = 0.01)
  expect_message(fit_wsgn(th), "one-sided")
  expect_error(fit_wsgn(th, se = "information"), "bound")
  ## the fit itself is the maximum on the boundary and is kept
  f <- suppressMessages(fit_wsgn(th))
  expect_gt(f$logLik, sum(dwsgn(th, 1, 3, 0.4, sigma = 0.01, log = TRUE)))
  ## an ordinary interior fit is not flagged
  set.seed(5)
  g <- fit_wsgn(rwsgn(100, 1, 1.8, 0.3, sigma = 0.5))
  expect_false(g$at_one_sided)
  expect_false(anyNA(g$se))
})

test_that("profile_wsgn takes its rate bound from beta_floor, as fit_wsgn does", {
  th <- 1 + (-3:3) * 1e-14
  p01 <- profile_wsgn(th, 2, beta_floor = 0.1)$logLik
  p1 <- profile_wsgn(th, 2, beta_floor = 1)$logLik
  ## floor 0.1 allows the bound 99, floor 1 only 30, whatever the profiled shape
  expect_gt(p01, 211)
  expect_gt(p01, p1 + 5)
  expect_equal(p1, profile_wsgn(th, 2, beta_floor = 1.5)$logLik, tolerance = 1e-9)
})

test_that("a resolution giving more than 1e6 cells is refused, in R and in C", {
  expect_error(grouped_loglik(c(1, 2, 3), 2 * pi / 3e9, 1, 0.4, 0.3, 1.5), "1e6")
  expect_error(fit_wsgn(c(1, 2, 3, 4), resolution = 2 * pi / 3e9), "1e6")
  expect_error(profile_wsgn(c(1, 2, 3, 4), 1, resolution = 2 * pi / 2e6), "1e6")
  expect_error(grouped_loglik(c(1, 2, 3), 2 * pi / 1.5e6, 1, 0.4, 0.3, 1.5), "1e6")
  ## the C entry points refuse it as well (NaN), without a cast out of range or an allocation
  out <- .C(wsgnstat:::C_wsgn_grouped_ll, theta = c(1, 2, 3), n = 3L, res = 2 * pi / 3e9,
            origin = 0, mu = 1, sigma = 0.4, beta = 1.5, lam = 0.3, K = 50L, out = double(1))$out
  expect_true(is.nan(out))
  r <- .C(wsgnstat:::C_wsgn_fit, theta = c(1, 2, 3, 4), n = 4L, grouped = 1L, res = 2 * pi / 3e9,
          origin = 0, beta0 = 0.001, sym = 0L, fixbeta = 0, window = 0.4, K = 50L, nthreads = 1L,
          out = double(12))$out
  expect_true(all(is.nan(r)))
  ## 1e6 cells itself is allowed
  expect_true(is.finite(grouped_loglik(c(1, 2, 3), 2 * pi / 1e6, 1, 0.4, 0.3, 1.5)))
})

test_that("an observation within rounding of the location is on the side it lies on", {
  ## 1 - 2^-53 and 1 + 2^-52 are the doubles next to 1; with rates near 5e19 the density a rounding
  ## step from the location is zero on either side, never the value at the location
  sg <- 1e-20
  for (x in c(1 - 2^-53, 1 + 2^-52))
    expect_lt(dwsgn(x, 1, 2, 0.5, sigma = sg, log = TRUE), -1e6)
  expect_gt(dwsgn(1, 1, 2, 0.5, sigma = sg, log = TRUE), 40)
  th <- 1 + (-3:3) * 1e-14
  expect_equal(profile_wsgn(th, 2, beta_floor = 0.1)$logLik, 211.601, tolerance = 1e-5)
  ## the profile at the fit's own shape is the fit
  f <- fit_wsgn(th, beta_floor = 0.1, se = "none")
  expect_equal(profile_wsgn(th, f$par[["beta"]], beta_floor = 0.1)$logLik, f$logLik,
               tolerance = 1e-9)
  expect_equal(f$logLik, 214.0974, tolerance = 1e-6)
})

test_that("a grouped one-sided fit is flagged, with no information errors and no bootstrap offer", {
  skip_on_cran()
  set.seed(204)
  th <- rwsgn(60, 1, 3, 0.5, sigma = 0.3)
  th <- 1 + abs(th - 1)
  h <- 5 * pi / 180
  expect_message(f <- fit_wsgn((round(th / h) * h) %% (2 * pi), resolution = h),
                 "one-sided")
  expect_true(f$at_one_sided)
  expect_gt(abs(diff(f$log_rates)), log(1e12))
  expect_null(f$se)
  expect_error(fit_wsgn((round(th / h) * h) %% (2 * pi), resolution = h, se = "information"),
               "one-sided")
  m <- tryCatch(fit_wsgn((round(th / h) * h) %% (2 * pi), resolution = h), message = function(m)
    conditionMessage(m))
  expect_false(grepl("bootstrap", m))
})

test_that("ngrid is capped, the lam error of a symmetric fit is NaN by design, at_upper is looser", {
  expect_error(fit_wsgn(c(1, 2, 3, 4), ngrid = 1e9), "1e8")
  expect_error(fit_wsgn(c(1, 2, 3, 4), ngrid = 2147483648), "1e8")
  skip_on_cran()
  set.seed(5)
  f <- fit_wsgn(rwsgn(80, 1, 1.8, 0.5, sigma = 0.5), symmetric = TRUE)
  expect_true(is.nan(f$se[["lam"]]))
  expect_false(anyNA(f$se[c("mu", "sigma", "beta")]))
  set.seed(101)
  th <- rwsgn(40, 1, 3, 0.5, sigma = 0.01)
  th <- 1 + abs(th - 1)
  g <- suppressMessages(fit_wsgn(th))
  expect_gt(g$par[["beta"]], 19.99)
  expect_true(g$at_upper)
})

test_that("a one-sided fit that stops at a cusp is polished to the maximum", {
  skip_on_cran()
  ## data to 15 significant digits: the final free polish used to collapse along mu, leaving a fit
  ## 0.0089 below the profile at the upper bound of the shape (150.3940 against 150.40293)
  set.seed(101)
  x <- rwsgn(40, 1, 3, 0.5, sigma = 0.01)
  th <- signif(1 + abs(x - 1), 15)
  f <- suppressMessages(fit_wsgn(th, se = "none"))
  expect_equal(f$logLik, 150.4029337676, tolerance = 1e-9)
  expect_true(f$at_upper)
  expect_true(f$at_one_sided)
  expect_equal(f$par[["beta"]], 20, tolerance = 1e-6)
  expect_gte(f$logLik, profile_wsgn(th, 20)$logLik - 1e-8)
})
