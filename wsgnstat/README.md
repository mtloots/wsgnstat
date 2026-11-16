
# wsgnstat

<!-- badges: start -->

[![R-CMD-check](https://github.com/mtloots/wsgnstat/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mtloots/wsgnstat/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/wsgnstat)](https://CRAN.R-project.org/package=wsgnstat)
<!-- badges: end -->

The **two-piece exponential power distribution** on the line and its
**wrapped version on the circle**, treated as a uniform scale mixture
with generalised gamma mixing. The package name keeps the authors’
earlier one, the wrapped skew generalised normal distribution. It
provides density, distribution and quantile functions and random
generation, the characteristic function and the trigonometric moments
for every shape parameter, and maximum likelihood fitting of the
four-parameter wrapped distribution under a lower bound on the shape, or
by the grouped likelihood for rounded data. The computation runs on a
pure-C back-end shared with the Python package of the same name.

## Installation

From CRAN:

``` r
install.packages("wsgnstat")
```

The development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("mtloots/wsgnstat", subdir = "wsgnstat")
```

## Example

Simulate from a skew member of the four-parameter wrapped distribution,
fit it by maximum likelihood, and compare with the truth; at a fitted
shape above one the standard errors come from the expected Fisher
information:

``` r
library(wsgnstat)
set.seed(1)
theta <- rwsgn(200, mu = 1, beta = 1.6, lam = 0.3, sigma = 0.4)
fit <- fit_wsgn(theta)
round(rbind(truth = c(1, 0.4, 0.3, 1.6), estimate = coef(fit), se = fit$se), 3)
#>             mu sigma   lam  beta
#> truth    1.000 0.400 0.300 1.600
#> estimate 1.095 0.402 0.318 1.732
#> se       0.106 0.054 0.046 0.263
```

The characteristic function is computed for every shape, so it equals
one at zero also for small shapes, and at shape 2 with `lam = 0.5` it is
the normal characteristic function `exp(-sigma^2 t^2)`:

``` r
sapply(c(0.1, 0.3, 1, 2), function(b) Mod(cfsgn(0, 0, 0.4, b, 0.3, 0.7) - 1))
#> [1] 0 0 0 0
t <- 1:4
rbind(engine = Re(cfsgn(t, 0, 0.4, 2, 0.5, 0.5)), closed_form = exp(-0.16 * t^2))
#>                  [,1]      [,2]      [,3]       [,4]
#> engine      0.8521438 0.5272924 0.2369278 0.07730474
#> closed_form 0.8521438 0.5272924 0.2369278 0.07730474
```

Rounded data are fitted by the grouped likelihood over the recording
cells; here the directions are recorded to 10 degrees:

``` r
h <- pi / 18
theta10 <- (round(theta / h) * h) %% (2 * pi)
coef(fit_wsgn(theta10, resolution = h, se = "none"))
#>        mu     sigma       lam      beta 
#> 1.0800109 0.3963098 0.3124758 1.7239484
```

The vignette, `vignette("wsgnstat")`, treats the mixture representation,
the characteristic function for small and large shapes, the accuracy of
the wrapped density, the fit with its shape floor on wind directions,
and the grouped likelihood for rounded ant directions.
