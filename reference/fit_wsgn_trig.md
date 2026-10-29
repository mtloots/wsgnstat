# Fit the wrapped skew generalised normal by trigonometric moments (C back-end)

Minimum-distance matching of the first two sample trigonometric moments
to their population counterparts, evaluated in C through the integral
characteristic function, hence available for every `beta > 0`. Standard
errors come from the minimum-distance sandwich.

## Usage

``` r
fit_wsgn_trig(theta, starts = NULL)
```

## Arguments

- theta:

  numeric vector of angles in radians.

- starts:

  list of numeric starting values `c(mu, beta, lam)`.

## Value

a list with the estimate `par = c(mu, beta, lam)`, its standard errors
`se`, and the objective value `value`.

## Examples

``` r
## simulate from a known member and recover it
set.seed(11)
truth <- c(mu = 1, beta = 1.5, lam = 0.3)
th <- rwsgn(1000, truth["mu"], truth["beta"], truth["lam"])
fit <- fit_wsgn_trig(th)
round(rbind(truth = truth, estimate = fit$par, se = fit$se), 3)
#>             mu  beta   lam
#> truth    1.000 1.500 0.300
#> estimate 1.005 1.496 0.288
#> se       0.079 0.081 0.030
```
