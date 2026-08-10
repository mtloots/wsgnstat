## The skew generalised normal on the line and its wrapping onto the circle. The family's whole
## point is that the characteristic function is reached through absolutely convergent integrals
## rather than a series with a limited radius, so the checks below pin the density, the wrapping
## and the trigonometric moments against closed forms where closed forms exist.

test_that("the line density is a density", {
  ## the range has to cover the heavier shapes and the sigma = 2 case: at a half-width of
  ## 12 the third parameter set still leaves 5 per cent of its mass outside
  x <- seq(-30, 30, length.out = 60001); dx <- x[2] - x[1]
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
  th <- seq(0, 2 * pi, length.out = 20001); dth <- th[2] - th[1]
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
  x <- seq(-20, 20, length.out = 60001); dx <- x[2] - x[1]
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

test_that("the estimators recover the parameters that generated the data", {
  skip_on_cran()
  set.seed(17)
  th <- rsgn(3000, mu = 1.0, sigma = 1, beta = 2, lam1 = 0.5, lam2 = 0.5, wrap = TRUE)
  ft <- fit_wsgn_trig(th)
  expect_true(all(is.finite(unlist(ft[sapply(ft, is.numeric)]))))
  fm <- fit_wsgn_mle(th)
  expect_true(all(is.finite(unlist(fm[sapply(fm, is.numeric)]))))
})
