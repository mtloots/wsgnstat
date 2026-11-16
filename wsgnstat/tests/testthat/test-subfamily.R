## The three-parameter reciprocal-rate subfamily fitters of version 0.1.0, kept for backward
## compatibility.

test_that("the subfamily estimators recover the parameters that generated the data", {
  truths <- list(c(mu = 1, beta = 1.5, lam = 0.3), c(mu = 4, beta = 0.8, lam = 0.6))
  for (tr in truths) {
    set.seed(17)
    th <- rwsgn(2000, tr["mu"], tr["beta"], tr["lam"])
    ft <- fit_wsgn_trig(th)
    expect_named(ft$par, c("mu", "beta", "lam"))
    expect_true(all(ft$se > 0))
    expect_true(all(abs(ft$par - tr) < 4 * ft$se))
    fm <- fit_wsgn_mle(th)
    expect_lt(abs(fm$par["mu"] - tr["mu"]), 0.1)
    expect_lt(abs(fm$par["beta"] - tr["beta"]) / tr["beta"], 0.15)
    expect_lt(abs(fm$par["lam"] - tr["lam"]), 0.06)
    nll_trig <- -sum(log(dwsgn(th, ft$par["mu"], ft$par["beta"], ft$par["lam"])))
    expect_gte(fm$logLik, -nll_trig - 1e-8)
  }
})

test_that("bootstrap standard errors for the subfamily fit are sensible", {
  set.seed(31)
  th <- rwsgn(400, 1, 1.5, 0.3)
  fm <- fit_wsgn_mle(th, B = 5)
  expect_named(fm$se, c("mu", "beta", "lam"))
  expect_true(all(is.finite(fm$se) & fm$se > 0))
})
