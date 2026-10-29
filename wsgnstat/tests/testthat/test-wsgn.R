## The skew generalised normal on the line and its wrapping onto the circle. The family's whole
## point is that the characteristic function is reached through absolutely convergent integrals
## rather than a series with a limited radius, so the checks below pin the density, the wrapping
## and the trigonometric moments against closed forms where closed forms exist.

test_that("the line density is a density", {
  ## the range has to cover the heavier shapes and the sigma = 2 case: at a half-width of
  ## 12 the third parameter set still leaves 5 per cent of its mass outside
  x <- seq(-30, 30, length.out = 60001)
  dx <- x[2] - x[1]
  for (par in list(c(0, 1, 2, 0.5, 0.5), c(0, 1, 1.3, 0.8, 0.3), c(1, 2, 3, 0.2, 0.9))) {
    f <- dsgn(x, par[1], par[2], par[3], par[4], par[5])
    expect_true(all(f >= 0))
    expect_lt(abs(sum(f) * dx - 1), 1e-4)
  }
})

test_that("equal skew weights give a density symmetric about its location", {
  x <- seq(-8, 8, length.out = 4001)
  expect_lt(max(abs(dsgn(x, 0, 1, 2, 0.5, 0.5) - dsgn(-x, 0, 1, 2, 0.5, 0.5))), 1e-12)
  ## and unequal weights genuinely break that symmetry
  expect_gt(max(abs(dsgn(x, 0, 1, 2, 0.8, 0.2) - dsgn(-x, 0, 1, 2, 0.8, 0.2))), 1e-3)
})

test_that("the circular density integrates to one over the circle", {
  th <- seq(0, 2 * pi, length.out = 20001)
  dth <- th[2] - th[1]
  for (par in list(c(0, 2, 0.5), c(1, 1.5, 0.3), c(0, 3, 0.7))) {
    d <- dwsgn(th, par[1], par[2], par[3])
    expect_true(all(d >= 0))
    expect_lt(abs(sum(d) * dth - 1), 1e-3)
  }
})

test_that("at beta = 2 the circular law IS the wrapped normal", {
  ## CLOSED-FORM ANCHOR. With beta = 2 and equal weights the generalised normal is the normal
  ## with sigma = 1/sqrt(2), so the wrapped family must reproduce the wrapped normal exactly.
  ## This checks the wrapping machinery against mathematics rather than against itself.
  th <- c(0.5, 2.0, 4.5, 6.0)
  wn <- vapply(th, function(t) sum(dnorm(t + 2 * pi * (-80:80), 0, 1 / sqrt(2))), 0)
  expect_lt(max(abs(dwsgn(th, 0, 2, 0.5) - wn)), 1e-9)
})

test_that("the trigonometric moments match their closed form at beta = 2", {
  ## for the wrapped normal with sigma = 1/sqrt(2) the p-th trigonometric moment is
  ## exp(-p^2 sigma^2 / 2) = exp(-p^2 / 4), real and positive because the law is symmetric.
  tm <- trig_moments(4, 0, 2, 0.5)
  expect_lt(max(abs(tm$A - exp(-(1:4)^2 / 4))), 1e-9)
  expect_lt(max(abs(tm$B)), 1e-9)
})

test_that("trigonometric moments are bounded by one in modulus", {
  for (par in list(c(0, 2, 0.5), c(0.7, 1.4, 0.25), c(0, 4, 0.9))) {
    tm <- trig_moments(3, par[1], par[2], par[3])
    expect_true(all(sqrt(tm$A^2 + tm$B^2) <= 1 + 1e-8))
  }
})

test_that("the characteristic function behaves like one", {
  expect_equal(cfsgn(0, 0, 1, 2, 0.5, 0.5), as.complex(1), tolerance = 1e-10)
  tt <- c(0.2, 0.7, 1.5, 3)
  for (par in list(c(0, 1, 2, 0.5, 0.5), c(0, 1, 1.3, 0.8, 0.3))) {
    cf <- cfsgn(tt, par[1], par[2], par[3], par[4], par[5])
    expect_true(all(Mod(cf) <= 1 + 1e-8))
    ## conjugate symmetry: phi(-t) = Conj(phi(t)), which any real random variable must satisfy
    cfm <- cfsgn(-tt, par[1], par[2], par[3], par[4], par[5])
    expect_lt(max(Mod(cfm - Conj(cf))), 1e-8)
  }
})

test_that("the characteristic function is the Fourier transform of the density", {
  ## an independent route to the same object: integrate exp(i t x) f(x) numerically and compare
  ## truncation, not resolution, sets the floor here: at a half-width of 20 the error is
  ## 1e-5 and does not move when the grid is refined sixteenfold
  x <- seq(-20, 20, length.out = 60001)
  dx <- x[2] - x[1]
  f <- dsgn(x, 0, 1, 1.3, 0.8, 0.3)
  for (t in c(0.4, 1.1)) {
    num <- sum(exp(1i * t * x) * f) * dx
    expect_lt(Mod(num - cfsgn(t, 0, 1, 1.3, 0.8, 0.3)), 1e-4)
  }
})

test_that("simulation follows the density", {
  set.seed(5)
  x <- rsgn(40000, 0, 1, 2, 0.5, 0.5)
  expect_true(all(is.finite(x)))
  expect_lt(abs(mean(x)), 0.05)              # symmetric case has mean zero
  ## wrapped draws land on the circle
  th <- rsgn(2000, 0, 1, 2, 0.5, 0.5, wrap = TRUE)
  expect_true(all(th >= 0 & th < 2 * pi + 1e-9))
})

test_that("the line distribution function is the integral of the density", {
  for (par in list(c(0, 1, 1.3, 0.3, 0.8), c(1, 2, 0.6, 0.8, 0.2), c(0, 1, 4, 0.5, 0.5))) {
    q <- c(-3, -0.5, par[1], 0.4, 2.5)
    num <- vapply(q, function(u) {
      stats::integrate(dsgn, -Inf, u, mu = par[1], sigma = par[2], beta = par[3],
                       lam1 = par[4], lam2 = par[5], rel.tol = 1e-11, subdivisions = 2000L)$value
    }, 0)
    expect_lt(max(abs(psgn(q, par[1], par[2], par[3], par[4], par[5]) - num)), 1e-9)
  }
  ## closed forms: the mass left of the location is lam1 / (lam1 + lam2), and the normal member
  expect_equal(psgn(0, 0, 1, 1.3, 0.3, 0.8), 0.3 / 1.1, tolerance = 1e-14)
  q <- c(-1.2, 0.4, 2)
  expect_lt(max(abs(psgn(q, 0, 0.5, 2, 0.5, 0.5) - pnorm(q, 0, sqrt(1 / 2)))), 1e-14)
  expect_identical(psgn(c(-Inf, Inf), 0, 1, 1.3, 0.3, 0.8), c(0, 1))
})

test_that("the wrapped distribution function is the integral of the wrapped density", {
  for (par in list(c(1, 1.5, 0.3), c(5, 0.8, 0.2), c(3, 8, 0.9))) {
    q <- c(0.3, 1, 2.5, 4, 6)
    num <- vapply(q, function(u) {
      stats::integrate(dwsgn, 0, u, mu = par[1], beta = par[2], lam = par[3],
                       rel.tol = 1e-11, subdivisions = 2000L)$value
    }, 0)
    expect_lt(max(abs(pwsgn(q, par[1], par[2], par[3]) - num)), 1e-9)
    expect_equal(pwsgn(c(-1, 0, 2 * pi, 7), par[1], par[2], par[3]), c(0, 0, 1, 1),
                 tolerance = 1e-12)
  }
  ## closed-form anchor: the wrapped normal with sigma = 1/sqrt(2)
  wn <- function(t) {
    sum(pnorm(t + 2 * pi * (-60:60), 0, sqrt(1 / 2)) - pnorm(2 * pi * (-60:60), 0, sqrt(1 / 2)))
  }
  expect_lt(max(abs(pwsgn(c(1, 3, 5), 0, 2, 0.5) - vapply(c(1, 3, 5), wn, 0))), 1e-12)
})

test_that("the quantile function inverts the distribution function", {
  p <- c(0.001, 0.05, 0.3, 0.5, 0.9, 0.999)
  for (par in list(c(1, 1.5, 0.3), c(5, 0.8, 0.2), c(0, 2, 0.5))) {
    q <- qwsgn(p, par[1], par[2], par[3])
    expect_true(all(q >= 0 & q <= 2 * pi))
    expect_lt(max(abs(pwsgn(q, par[1], par[2], par[3]) - p)), 1e-12)
  }
  expect_identical(qwsgn(c(0, 1), 1, 1.5, 0.3), c(0, 2 * pi))
})

test_that("a truncated wrapping sum is reported, not hidden", {
  ## at beta = 0.3 the default K = 50 misses mass; pwsgn(2 pi) says how much, raising K
  ## recovers it, and qwsgn warns rather than returning a meaningless angle
  tot <- pwsgn(2 * pi, 0, 0.3, 0.5)
  expect_lt(tot, 0.95)
  expect_gt(pwsgn(2 * pi, 0, 0.3, 0.5, K = 400), tot)
  expect_warning(qwsgn(0.99, 0, 0.3, 0.5), "raise K")
  expect_no_warning(qwsgn(0.5, 0, 0.3, 0.5))
})

test_that("wrapped draws follow the wrapped law", {
  set.seed(23)
  th <- rwsgn(20000, mu = 1, beta = 1.5, lam = 0.3)
  expect_true(all(th >= 0 & th < 2 * pi))
  ## Kolmogorov-Smirnov distance from 0 against the exact wrapped distribution function
  ## (the 0.1 per cent critical value at n = 20000 is about 1.95 / sqrt(n) = 0.0138)
  u <- sort(pwsgn(th, mu = 1, beta = 1.5, lam = 0.3))
  n <- length(u)
  ks <- max(pmax((1:n) / n - u, u - (0:(n - 1)) / n))
  expect_lt(ks, 1.95 / sqrt(n))
  ## first trigonometric moment against the population value
  tm <- trig_moments(1, 1, 1.5, 0.3)
  expect_lt(abs(mean(cos(th)) - tm$A), 4 * sd(cos(th)) / sqrt(n))
  expect_lt(abs(mean(sin(th)) - tm$B), 4 * sd(sin(th)) / sqrt(n))
  ## and the line sampler puts lam1 / (lam1 + lam2) of its mass left of the location
  x <- rsgn(20000, 0, 1, 1.3, 0.3, 0.8)
  expect_lt(abs(mean(x < 0) - 0.3 / 1.1), 4 * sqrt(0.25 / 20000))
})

test_that("the estimators recover the parameters that generated the data", {
  ## Draws from known members of the identified family; the moment estimator must land within
  ## four of its own sandwich standard errors and the likelihood estimator close to the truth.
  truths <- list(c(mu = 1, beta = 1.5, lam = 0.3), c(mu = 4, beta = 0.8, lam = 0.6),
                 c(mu = 2, beta = 3, lam = 0.5))
  for (tr in truths) {
    set.seed(17)
    th <- rwsgn(2000, tr["mu"], tr["beta"], tr["lam"])
    ft <- fit_wsgn_trig(th)
    expect_named(ft$par, c("mu", "beta", "lam"))
    expect_true(all(ft$se > 0))
    expect_true(all(abs(ft$par - tr) < 4 * ft$se))
    fm <- fit_wsgn_mle(th)
    expect_lt(abs(fm$par["mu"] - tr["mu"]), 0.1)
    expect_lt(abs(fm$par["beta"] - tr["beta"]) / tr["beta"], 0.15)  # MC sd 3.5% at beta 3
    expect_lt(abs(fm$par["lam"] - tr["lam"]), 0.06)
    ## the maximum likelihood fit is at least as likely as the moment fit it starts from
    nll_trig <- -sum(log(dwsgn(th, ft$par["mu"], ft$par["beta"], ft$par["lam"])))
    expect_gte(fm$logLik, -nll_trig - 1e-8)
  }
})

test_that("bootstrap standard errors for the likelihood fit are sensible", {
  skip_on_cran()
  set.seed(31)
  th <- rwsgn(400, 1, 1.5, 0.3)
  fm <- fit_wsgn_mle(th, B = 20)
  expect_named(fm$se, c("mu", "beta", "lam"))
  expect_true(all(is.finite(fm$se) & fm$se > 0))
  ## of the same order as the moment estimator's sandwich errors
  ft <- fit_wsgn_trig(th)
  expect_true(all(fm$se < 5 * ft$se & fm$se > ft$se / 5))
})
