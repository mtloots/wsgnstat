## The two-piece exponential power distribution on the line, and its uniform scale mixture
## representation: X = mu + Z V, Z^beta ~ Gamma(1 + 1/beta), V uniform on [-sigma/lam2, sigma/lam1].

test_that("the line density is a density, with the normal member as a closed-form anchor", {
  for (par in list(c(0, 1, 2, 0.5, 0.5), c(0, 1, 1.3, 0.8, 0.3), c(1, 2, 0.4, 0.2, 0.9))) {
    tot <- stats::integrate(dsgn, -Inf, Inf, mu = par[1], sigma = par[2], beta = par[3],
                            lam1 = par[4], lam2 = par[5], rel.tol = 1e-10, subdivisions = 2000L)
    expect_equal(tot$value, 1, tolerance = 1e-8)
  }
  x <- seq(-4, 4, length.out = 81)
  expect_lt(max(abs(dsgn(x, 0, 0.5, 2, 0.5, 0.5) - dnorm(x, 0, sqrt(1 / 2)))), 1e-15)
  ## the asymmetric Laplace member, beta = 1
  expect_equal(dsgn(c(-1, 2), 0, 0.4, 1, 0.3, 0.7),
               0.21 / 0.4 * exp(-c(0.7 * 1, 0.3 * 2) / 0.4), tolerance = 1e-14)
  expect_lt(max(abs(dsgn(x, 0, 1, 2, 0.5, 0.5) - dsgn(-x, 0, 1, 2, 0.5, 0.5))), 1e-15)
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
  expect_equal(psgn(0, 0, 1, 1.3, 0.3, 0.8), 0.3 / 1.1, tolerance = 1e-14)
  q <- c(-1.2, 0.4, 2)
  expect_lt(max(abs(psgn(q, 0, 0.5, 2, 0.5, 0.5) - pnorm(q, 0, sqrt(1 / 2)))), 1e-14)
  expect_identical(psgn(c(-Inf, Inf), 0, 1, 1.3, 0.3, 0.8), c(0, 1))
  expect_true(is.na(psgn(NA, 0, 1, 1.3, 0.3, 0.8)))
})

test_that("the mixture representation reproduces the distribution, its moments and P(X < mu)", {
  ## draws from the single-uniform mixture pass a Kolmogorov-Smirnov test against the closed-form
  ## distribution function, at a small shape (cusp) and at a large one
  set.seed(101)
  for (par in list(c(0.3, 0.4, 0.3), c(2.5, 1, 0.7))) {
    be <- par[1]
    s <- par[2]
    l <- par[3]
    x <- rsgn(20000, 0, s, be, l, 1 - l)
    u <- sort(psgn(x, 0, s, be, l, 1 - l))
    n <- length(u)
    ks <- max(pmax((1:n) / n - u, u - (0:(n - 1)) / n))
    expect_lt(ks, 1.95 / sqrt(n))
    ## the mass to the left of the location is lam
    expect_lt(abs(mean(x < 0) - l), 4 * sqrt(l * (1 - l) / n))
  }
  ## the moments E(X - mu)^k = E Z^k E V^k, with a = lam / sigma and b = (1 - lam) / sigma
  be <- 1.5
  s <- 0.6
  l <- 0.3
  a <- l / s
  b <- (1 - l) / s
  ss <- 1 / a + 1 / b
  mom <- function(k) {
    gamma(1 + (k + 1) / be) / gamma(1 + 1 / be) * (a^-(k + 1) - (-b)^-(k + 1)) / ((k + 1) * ss)
  }
  for (k in 1:3) {
    num <- stats::integrate(function(x) x^k * dsgn(x, 0, s, be, l, 1 - l), -Inf, Inf,
                            rel.tol = 1e-11)$value
    expect_equal(num, mom(k), tolerance = 1e-8)
  }
  set.seed(102)
  x <- rsgn(1e5, 0, s, be, l, 1 - l)
  expect_lt(abs(mean(x) - mom(1)), 4 * sd(x) / sqrt(length(x)))
  expect_lt(abs(mean(x^2) - mom(2)), 4 * sd(x^2) / sqrt(length(x)))
  ## wrapped draws land on the circle
  expect_true(all(rsgn(1000, 0, 1, 0.5, 0.5, 0.5, wrap = TRUE) < 2 * pi))
})
