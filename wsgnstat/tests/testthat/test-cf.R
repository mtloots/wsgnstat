## The characteristic function for every shape: against 120-digit reference values from the
## convergent series (inst/extdata/ref_cf.csv, made by tools/gen_ref_cf.py in the repository),
## closed forms at shapes one and two, the value one at zero, and the decay constant of the
## trigonometric moments.

ref <- utils::read.csv(system.file("extdata", "ref_cf.csv", package = "wsgnstat"))

test_that("the characteristic function equals one at zero for every shape", {
  for (be in c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1, 1.5, 2, 3, 10)) {
    expect_lt(Mod(cfsgn(0, 0, 0.4, be, 0.3, 0.7) - 1), 1e-14)
  }
  ## the reciprocal-rate member of version 0.1.0 at small shapes as well
  for (be in c(0.1, 0.2, 0.3)) expect_lt(Mod(cfsgn(0, 0, sqrt(0.21), be, 0.3, 0.7) - 1), 1e-14)
})

test_that("the characteristic function matches the 120-digit series values for shapes 0.1 to 3", {
  sig <- ifelse(ref$case == 1, 0.4, sqrt(0.21))
  z <- mapply(function(t, b, s) cfsgn(t, 0, s, b, 0.3, 0.7), ref$t, ref$beta, sig)
  zr <- complex(real = ref$re, imaginary = ref$im)
  expect_lt(max(Mod(z - zr) / Mod(zr)), 1e-12)
  ## the trigonometric moments of the reciprocal-rate member are the case-2 values at p = 1, 2
  r2 <- ref[ref$case == 2 & ref$t <= 2 & ref$beta == 0.1, ]
  tm <- trig_moments(2, 0, 0.1, 0.3)
  zr2 <- complex(real = r2$re, imaginary = r2$im)
  expect_lt(max(Mod(complex(real = tm$A, imaginary = tm$B) - zr2)), 1e-20)
})

test_that("the series method agrees where it converges, and refuses where it cancels", {
  for (be in c(0.1, 0.3, 0.6, 1, 2, 3)) {
    t <- c(0.5, 1, 2, 3)
    zs <- suppressWarnings(cfsgn(t, 0.2, 0.4, be, 0.3, 0.7, method = "series"))
    ok <- !is.na(Re(zs))
    expect_gte(sum(ok), 2)
    expect_lt(max(Mod(zs[ok] - cfsgn(t[ok], 0.2, 0.4, be, 0.3, 0.7))), 1e-12)
  }
  expect_warning(z <- cfsgn(20, 0, 0.4, 1.2, 0.3, 0.7, method = "series"), "cancels")
  expect_true(is.na(Re(z)))
  ## a sum stopped inside a cancellation once returned about 5e215 here
  expect_warning(z <- cfsgn(5, 0, 0.6, 1.5, 0.8, 0.2, method = "series"), "cancels")
  expect_true(is.na(Re(z)))
  ## at small shapes the series must not be truncated early (the value is of order 1e-20 there)
  t <- c(1, 2, 3)
  zs <- cfsgn(t, 0, 0.4, 0.05, 0.3, 0.7, method = "series")
  zq <- cfsgn(t, 0, 0.4, 0.05, 0.3, 0.7)
  expect_lt(max(Mod(zs - zq) / Mod(zq)), 1e-11)
  ## one half only (lam2 huge), where the ascending series cancels beyond tol: refused
  expect_true(is.na(Re(suppressWarnings(cfsgn(10, 0, 1, 10, 1, 1e300, method = "series")))))
  ## over a grid of shapes, skewness, scales and arguments, the series either agrees with the
  ## quadrature to 1e-8 or is refused
  grid <- expand.grid(be = c(0.1, 0.3, 0.7, 1, 1.1, 1.5, 2, 3, 5), l = c(0.1, 0.5, 0.8, 0.95),
                      s = c(0.05, 0.6, 2))
  dev <- mapply(function(be, l, s) {
    tt <- c(0.5, 2, 5, 30)
    zs <- suppressWarnings(cfsgn(tt, 0.3, s, be, l, 1 - l, method = "series"))
    zq <- cfsgn(tt, 0.3, s, be, l, 1 - l)
    ok <- !is.na(Re(zs))
    c(max(c(0, Mod(zs[ok] - zq[ok]))), sum(ok))
  }, grid$be, grid$l, grid$s)
  expect_lt(max(dev[1, ]), 1e-8)
  expect_gt(sum(dev[2, ]), 0.5 * 4 * nrow(grid))
  expect_error(cfsgn(1, tol = -1), "tol")
  expect_error(cfsgn(1, tol = 0), "tol")
  expect_error(cfsgn(1, tol = Inf), "tol")
  expect_error(trig_moments(2, tol = 0), "tol")
  expect_error(trig_moments(NA), "P")
  expect_warning(z <- cfsgn(0.02, 0, 0.5, 50, 0.5, 0.5, tol = 1e-16), "1e-14")
  ## against the 60-digit ascending series, 0.99993468977084714 (w = 0.02, beta = 50)
  expect_lt(abs(Re(z) - 0.99993468977084714), 5e-14)
})

test_that("closed forms: the asymmetric Laplace (beta = 1) and the normal (beta = 2) members", {
  t <- c(0.3, 1, 2.5, 7, 40)
  a <- 0.3 / 0.4
  b <- 0.7 / 0.4
  expect_lt(max(Mod(cfsgn(t, 0.5, 0.4, 1, 0.3, 0.7) -
                      exp(0.5i * t) * a * b / ((a - 1i * t) * (b + 1i * t)))), 1e-14)
  expect_lt(max(Mod(cfsgn(t, 0, 0.4, 2, 0.5, 0.5) - exp(-0.16 * t^2))), 1e-14)
  ## conjugate symmetry phi(-t) = Conj(phi(t)) and |phi| <= 1
  z <- cfsgn(t, 0.3, 0.4, 0.7, 0.3, 0.7)
  expect_lt(max(Mod(cfsgn(-t, 0.3, 0.4, 0.7, 0.3, 0.7) - Conj(z))), 1e-15)
  expect_true(all(Mod(z) <= 1))
})

test_that("the characteristic function is the Fourier transform of the density", {
  ## integrated piece by piece over [-400, 400], beyond which the mass is below 1e-15
  edges <- seq(-400, 400, by = 5)
  for (par in list(c(1.3, 0.8, 0.3), c(0.7, 0.6, 0.5))) {
    for (t in c(0.4, 1.1)) {
      part <- function(fun) {
        sum(vapply(seq_len(length(edges) - 1L), function(i) {
          stats::integrate(function(x) fun(exp(1i * t * x) * dsgn(x, 0, 1, par[1], par[2], par[3])),
                           edges[i], edges[i + 1L], rel.tol = 1e-12, subdivisions = 200L)$value
        }, 0))
      }
      num <- complex(real = part(Re), imaginary = part(Im))
      expect_lt(Mod(num - cfsgn(t, 0, 1, par[1], par[2], par[3])), 1e-10)
    }
  }
})

test_that("trigonometric moments: wrapped normal anchor, simulation and the decay constant", {
  tm <- trig_moments(5, 0, 2, 0.5, sigma = 0.6)
  expect_lt(max(abs(tm$A - exp(-0.36 * (1:5)^2))), 1e-14)
  expect_lt(max(abs(tm$B)), 1e-15)
  ## default sigma: the reciprocal-rate subfamily of version 0.1.0
  expect_identical(trig_moments(2, 1, 1.5, 0.3), trig_moments(2, 1, 1.5, 0.3, sigma = sqrt(0.21)))
  ## |phi_p| p^(1 + beta) / C tends to one; the relative remainder is O(p^-beta)
  for (be in c(0.8, 1.5, 3)) {
    dc <- trig_decay(be, 0.3, 0.4)
    expect_identical(unname(dc["exponent"]), 1 + be)
    tm <- trig_moments(500, 0, be, 0.3, sigma = 0.4)
    p <- c(50, 500)
    r <- Mod(complex(real = tm$A[p], imaginary = tm$B[p])) * p^(1 + be) / dc["C"]
    expect_lt(abs(r[2] - 1), 0.02)
    expect_lt(abs(r[2] - 1), abs(r[1] - 1))
  }
  expect_identical(unname(trig_decay(2, 0.5, 0.4)["C"]), 0)
  expect_identical(unname(trig_decay(4, 0.5, 0.4)["C"]), 0)
  ## the first moment against a large simulated sample
  set.seed(7)
  th <- rwsgn(1e5, mu = 1, beta = 0.6, lam = 0.3, sigma = 0.4)
  tm1 <- trig_moments(1, mu = 1, beta = 0.6, lam = 0.3, sigma = 0.4)
  expect_lt(abs(mean(cos(th)) - tm1$A), 4 * sd(cos(th)) / sqrt(1e5))
  expect_lt(abs(mean(sin(th)) - tm1$B), 4 * sd(sin(th)) / sqrt(1e5))
})
