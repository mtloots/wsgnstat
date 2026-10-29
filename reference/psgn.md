# Distribution function of the skew generalised normal distribution on the line (C back-end)

Evaluated in closed form through the regularised incomplete gamma
function: the mass to the left of `mu` is `lam1 / (lam1 + lam2)`, and
each tail is a scaled gamma tail in `(lam * |x - mu| / sigma)^beta`.

## Usage

``` r
psgn(q, mu = 0, sigma = 1, beta = 2, lam1 = 0.5, lam2 = 0.5)
```

## Arguments

- q:

  numeric vector of quantiles.

- mu:

  location parameter.

- sigma:

  scale parameter (\> 0).

- beta:

  shape (tail) parameter (\> 0).

- lam1:

  right-hand rate parameter (\> 0).

- lam2:

  left-hand rate parameter (\> 0).

## Value

the distribution function evaluated at `q`.

## Examples

``` r
## agreement with numerical integration of the density
q <- c(-1, 0, 0.7, 2.5)
closed <- psgn(q, 0, 1, 1.3, 0.3, 0.8)
numeric <- vapply(q, function(u)
  integrate(dsgn, -Inf, u, mu = 0, sigma = 1, beta = 1.3, lam1 = 0.3, lam2 = 0.8)$value, 0)
rbind(closed, numeric)
#>              [,1]      [,2]      [,3]      [,4]
#> closed  0.0978961 0.2727273 0.4290232 0.7197377
#> numeric 0.0978961 0.2727273 0.4290232 0.7197375
## the normal member: psgn equals pnorm with variance 1/2
psgn(0.4, 0, 0.5, 2, 0.5, 0.5) - pnorm(0.4, 0, sqrt(1 / 2))
#> [1] 1.110223e-16
```
