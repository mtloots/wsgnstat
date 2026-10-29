# Fit the wrapped skew generalised normal by maximum likelihood (C back-end)

Nelder-Mead over the identified parameters, with the wrapped
log-likelihood evaluated in C. Standard errors, when requested, come
from a parametric bootstrap, the observed information being unstable
near the interior cusp the density develops as `beta` approaches one.

## Usage

``` r
fit_wsgn_mle(theta, start = NULL, K = 50, B = 0)
```

## Arguments

- theta:

  numeric vector of angles in radians.

- start:

  numeric starting value `c(mu, beta, lam)`; default is the
  trigonometric-moment estimate.

- K:

  truncation of the wrapping sum.

- B:

  number of parametric-bootstrap replicates for standard errors; `0`
  skips them.

## Value

a list with `par`, `logLik`, and (if `B > 0`) `se`.

## Examples

``` r
set.seed(12)
truth <- c(mu = 1, beta = 1.5, lam = 0.3)
th <- rwsgn(500, truth["mu"], truth["beta"], truth["lam"])
fit <- fit_wsgn_mle(th)
round(rbind(truth = truth, estimate = fit$par), 3)
#>             mu  beta   lam
#> truth    1.000 1.500 0.300
#> estimate 1.048 1.499 0.303
fit$logLik
#> [1] -663.4528
# \donttest{
## bootstrap standard errors refit the model once per replicate (a few seconds for B = 30)
fit_b <- fit_wsgn_mle(th, B = 30)
round(rbind(estimate = fit_b$par, se = fit_b$se), 3)
#>             mu  beta   lam
#> estimate 1.048 1.499 0.303
#> se       0.079 0.106 0.029
# }
```
