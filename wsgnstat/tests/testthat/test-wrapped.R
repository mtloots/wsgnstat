## The wrapped density, distribution and quantile functions, against 40-digit image sums
## (inst/extdata/ref_dw.csv, made by tools/gen_ref_dw.py in the repository with mpmath's
## Euler-Maclaurin summation) and closed forms.

ref <- utils::read.csv(system.file("extdata", "ref_dw.csv", package = "wsgnstat"))

test_that("the wrapped density matches the 40-digit image sums for shapes 0.1 to 2", {
  sig <- ifelse(ref$case == 1, 0.4, sqrt(0.21))
  v <- mapply(function(th, b, s) dwsgn(th, 0.7, b, 0.3, sigma = s), ref$theta, ref$beta, sig)
  expect_lt(max(abs(v / ref$dens - 1)), 1e-12)
  ## the number of images summed directly does not matter: the tail is exact
  th <- seq(0, 2 * pi, length.out = 13)
  for (be in c(0.1, 0.3, 1)) {
    expect_lt(max(abs(dwsgn(th, 0.7, be, 0.3, K = 3, sigma = 0.4) /
                        dwsgn(th, 0.7, be, 0.3, K = 400, sigma = 0.4) - 1)), 1e-12)
  }
  ## log = TRUE, and the default sigma is the reciprocal-rate value
  expect_equal(dwsgn(th, 1, 0.5, 0.3, log = TRUE), log(dwsgn(th, 1, 0.5, 0.3, sigma = sqrt(0.21))),
               tolerance = 1e-14)
})

test_that("the wrapped density integrates to one, also at very small shapes", {
  for (be in c(0.05, 0.2, 1, 4)) {
    tot <- stats::integrate(dwsgn, 0, 2 * pi, mu = 1, beta = be, lam = 0.3, sigma = 0.4,
                            rel.tol = 1e-12, subdivisions = 2000L)$value
    expect_equal(tot, 1, tolerance = 1e-10)
  }
})

test_that("at beta = 2 and lam = 0.5 the wrapped distribution is the wrapped normal", {
  th <- c(0.5, 2.0, 4.5, 6.0)
  wn <- vapply(th, function(t) sum(dnorm(t + 2 * pi * (-80:80), 0, sqrt(2) * 0.6)), 0)
  expect_lt(max(abs(dwsgn(th, 0, 2, 0.5, sigma = 0.6) - wn)), 1e-15)
  wnp <- function(t) {
    sum(pnorm(t + 2 * pi * (-60:60), 0, sqrt(2) * 0.6) - pnorm(2 * pi * (-60:60), 0, sqrt(2) * 0.6))
  }
  expect_lt(max(abs(pwsgn(c(1, 3, 5), 0, 2, 0.5, sigma = 0.6) - vapply(c(1, 3, 5), wnp, 0))), 1e-13)
})

test_that("the wrapped distribution function integrates the density, and qwsgn inverts it", {
  for (par in list(c(1, 1.5, 0.3, 0.5), c(5, 0.2, 0.2, 0.4), c(3, 8, 0.9, 1))) {
    q <- c(0.3, 1, 2.5, 4, 6)
    num <- vapply(q, function(u) {
      stats::integrate(dwsgn, 0, u, mu = par[1], beta = par[2], lam = par[3], sigma = par[4],
                       rel.tol = 1e-12, subdivisions = 2000L)$value
    }, 0)
    expect_lt(max(abs(pwsgn(q, par[1], par[2], par[3], sigma = par[4]) - num)), 1e-10)
    expect_identical(pwsgn(c(-1, 0, 2 * pi, 7), par[1], par[2], par[3], sigma = par[4]),
                     c(0, 0, 1, 1))
    p <- c(0.001, 0.05, 0.3, 0.5, 0.9, 0.999)
    qq <- qwsgn(p, par[1], par[2], par[3], sigma = par[4])
    expect_true(all(qq >= 0 & qq <= 2 * pi))
    expect_lt(max(abs(pwsgn(qq, par[1], par[2], par[3], sigma = par[4]) - p)), 1e-12)
  }
  ## the full mass is reached for the default K at a very small shape (truncated in 0.1.0)
  expect_equal(pwsgn(2 * pi - 1e-9, 0, 0.3, 0.5), 1, tolerance = 1e-8)
  expect_identical(qwsgn(c(0, 1), 1, 1.5, 0.3), c(0, 2 * pi))
  expect_true(is.na(qwsgn(NA, 1, 1.5, 0.3)))
  expect_true(is.na(pwsgn(NA, 1, 1.5, 0.3)))
})

test_that("wrapped draws follow the wrapped distribution", {
  set.seed(23)
  for (par in list(c(1, 1.5, 0.3, 0.5), c(4, 0.4, 0.6, 0.3))) {
    th <- rwsgn(if (par[2] < 1) 3000 else 20000, par[1], par[2], par[3], sigma = par[4])
    expect_true(all(th >= 0 & th < 2 * pi))
    u <- sort(pwsgn(th, par[1], par[2], par[3], sigma = par[4]))
    n <- length(u)
    ks <- max(pmax((1:n) / n - u, u - (0:(n - 1)) / n))
    expect_lt(ks, 1.95 / sqrt(n))
  }
})

test_that("very small shapes and the lower tail keep their accuracy", {
  ## at beta = 0.005 and 0.001 the wrapped density is 1/(2 pi) to double precision
  for (be in c(0.005, 0.001)) {
    d <- dwsgn(c(0.5, 3), mu = 1, beta = be, lam = 0.3, sigma = 0.4)
    expect_lt(max(abs(d * 2 * pi - 1)), 1e-12)
  }
  ## small probabilities: p(q) / (q g(0)) -> 1 and the quantile function inverts it relatively
  g0 <- dwsgn(0, 1, 1.5, 0.3, sigma = 0.4)
  for (q in c(1e-8, 1e-12, 1e-20, 1e-200)) {
    p <- pwsgn(q, 1, 1.5, 0.3, sigma = 0.4)
    expect_lt(abs(p / (q * g0) - 1), 2 * q + 1e-13)
    expect_lt(abs(qwsgn(p, 1, 1.5, 0.3, sigma = 0.4) / q - 1), 1e-13)
  }
  expect_equal(qwsgn(1e-300, 1, 1.5, 0.3, sigma = 0.4), 1e-300 / g0, tolerance = 1e-12)
  expect_gt(qwsgn(1.6e-301, mu = 2, beta = 0.3, lam = 0.3, sigma = 0.4), 1e-301)
  ## next to the location at large shapes the argument (rate x)^beta would underflow: the arc
  ## [0, q] then has mass q g(0) to first order, g the wrapped density
  for (be in c(2, 5, 20, 50)) {
    q <- c(1e-200, 1e-100, 1e-20, 1e-12)
    p <- pwsgn(q, 0, be, 0.5, sigma = 1)
    expect_lt(max(abs(p / (q * dwsgn(0, 0, be, 0.5, sigma = 1)) - 1)), 1e-12)
    expect_lt(max(abs(qwsgn(p, 0, be, 0.5, sigma = 1) / q - 1)), 1e-12)
  }
  ## and the line distribution function is not flat next to the median at shape 50
  expect_gt(psgn(1 + 7e-9, 1, 0.01, 50, 0.5, 0.5) - 0.5, 1e-7)
  ## the upper end
  ## (the quantile near 2 pi is limited by the spacing of doubles there, 8.9e-16)
  expect_lt(abs(2 * pi - qwsgn(1 - 1e-14, 1, 1.5, 0.3, sigma = 0.4) - 1e-14 / g0), 2e-15)
})

test_that("invalid parameters are refused", {
  expect_error(dsgn(1, beta = -1), "positive")
  expect_error(cfsgn(1, lam1 = 0), "positive")
  expect_error(dwsgn(1, 0, -1, 0.5), "beta")
  expect_error(dwsgn(1, 0, 1, 1.5), "lam")
  expect_error(dwsgn(1, 0, 1, 0.5, sigma = 0), "sigma")
  expect_error(pwsgn(1, 0, 1, 0.5, K = 0), "K")
})
