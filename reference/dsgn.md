# Density of the skew generalised normal distribution on the line (C back-end)

The scale mixture of a uniform and a generalised gamma variable, with
two-piece exponential-power form: left rate `lam2`, right rate `lam1`.

## Usage

``` r
dsgn(x, mu = 0, sigma = 1, beta = 2, lam1 = 0.5, lam2 = 0.5)
```

## Arguments

- x:

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

the probability density function evaluated at `x`.

## Examples

``` r
## a right-skew, heavy-tailed member and the normal member it generalises
x <- seq(-4, 6, length.out = 401)
f_skew <- dsgn(x, mu = 0, sigma = 1, beta = 1.3, lam1 = 0.3, lam2 = 0.8)
f_norm <- dsgn(x, mu = 0, sigma = 0.5, beta = 2, lam1 = 0.5, lam2 = 0.5)
plot(x, f_skew, type = "l", ylab = "density")
lines(x, f_norm, lty = 2)

## at beta = 2 with lam1 = lam2 = 0.5 and sigma = 0.5 the law is N(0, 1/2)
max(abs(f_norm - dnorm(x, 0, sqrt(1 / 2))))
#> [1] 1.110223e-16
## the density integrates to one (Riemann sum on a wide grid)
xx <- seq(-40, 40, length.out = 80001)
sum(dsgn(xx, 0, 1, 1.3, 0.3, 0.8)) * (xx[2] - xx[1])
#> [1] 1
```
