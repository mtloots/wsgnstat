
# wsgnstat

[![R-CMD-check](https://github.com/mtloots/wsgnstat/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mtloots/wsgnstat/actions/workflows/R-CMD-check.yaml)

The **skew generalised normal** on the line and **wrapped onto the
circle**, with the characteristic function reached through absolutely
convergent integrals rather than a series.

Computation runs on a shared pure-C back-end that is also bound from
Python, and the two front ends are checked against each other value by
value.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("mtloots/wsgnstat", subdir = "wsgnstat")
```

## Usage

``` r
library(wsgnstat)
## at shape 2 with equal weights the circular law IS the wrapped normal, exactly
tm <- trig_moments(4, mu = 0, beta = 2, lam = 0.5)
rbind(engine = tm$A, closed_form = exp(-(1:4)^2 / 4))
#>                  [,1]      [,2]      [,3]       [,4]
#> engine      0.7788008 0.3678794 0.1053992 0.01831564
#> closed_form 0.7788008 0.3678794 0.1053992 0.01831564
```
