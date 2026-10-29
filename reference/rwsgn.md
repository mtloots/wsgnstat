# Random generation from the wrapped skew generalised normal distribution

Draws from the circular law in the identified parameterisation of
[`dwsgn`](https://mtloots.github.io/wsgnstat/reference/dwsgn.md), that
is
`rsgn(n, mu, sigma = sqrt(lam * (1 - lam)), beta, lam1 = lam, lam2 = 1 - lam, wrap = TRUE)`.
The estimators
[`fit_wsgn_trig`](https://mtloots.github.io/wsgnstat/reference/fit_wsgn_trig.md)
and
[`fit_wsgn_mle`](https://mtloots.github.io/wsgnstat/reference/fit_wsgn_mle.md)
estimate exactly these `(mu, beta, lam)`.

## Usage

``` r
rwsgn(n, mu = 0, beta = 2, lam = 0.5)
```

## Arguments

- n:

  number of draws.

- mu:

  location (junction) parameter; for a skew law this is not the mean
  direction.

- beta:

  shape (tail) parameter (\> 0).

- lam:

  skewness parameter in (0, 1); `0.5` is symmetric.

## Value

a numeric vector of angles in `[0, 2 * pi)`.

## Examples

``` r
set.seed(3)
th <- rwsgn(2000, mu = 1, beta = 1.5, lam = 0.3)
hist(th, breaks = seq(0, 2 * pi, length.out = 49), freq = FALSE, main = "",
     xlab = "angle (radians)")
curve(dwsgn(x, mu = 1, beta = 1.5, lam = 0.3), 0, 2 * pi, add = TRUE)
```
