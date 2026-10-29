# wsgnstat

The **skew generalised normal** distribution on the line and its
**wrapping onto the circle**: density, distribution and quantile
functions, random generation, characteristic function, trigonometric
moments, and estimation by trigonometric moments or maximum likelihood.
The characteristic function is reached through absolutely convergent
integrals rather than a series, so every shape parameter is available.
The computation runs on a pure-C back-end shared with the Python package
of the same name.

## Installation

From CRAN, once accepted:

``` r

install.packages("wsgnstat")
```

The development version from GitHub:

``` r

# install.packages("remotes")
remotes::install_github("mtloots/wsgnstat", subdir = "wsgnstat")
```

## Example

Simulate from a skew member of the wrapped family, estimate its
parameters, and compare with the truth:

``` r

library(wsgnstat)
set.seed(1)
truth <- c(mu = 1, beta = 1.5, lam = 0.3)
theta <- rwsgn(1000, truth["mu"], truth["beta"], truth["lam"])
fit <- fit_wsgn_mle(theta)
round(rbind(truth = truth, estimate = fit$par), 3)
#>          mu beta   lam
#> truth     1 1.50 0.300
#> estimate  1 1.48 0.315
```

At shape 2 with `lam = 0.5` the circular law is the wrapped normal with
sigma = 1/sqrt(2), and its trigonometric moments are exactly
exp(-p^2/4):

``` r

tm <- trig_moments(4, mu = 0, beta = 2, lam = 0.5)
rbind(engine = tm$A, closed_form = exp(-(1:4)^2 / 4))
#>                  [,1]      [,2]      [,3]       [,4]
#> engine      0.7788008 0.3678794 0.1053992 0.01831564
#> closed_form 0.7788008 0.3678794 0.1053992 0.01831564
```

The vignette,
[`vignette("wsgnstat")`](https://mtloots.github.io/wsgnstat/articles/wsgnstat.md),
walks through the family, its reductions, estimation and an application
to wind directions.
