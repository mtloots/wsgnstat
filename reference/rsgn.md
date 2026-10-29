# Random generation from the (wrapped) skew generalised normal distribution

Uses the stochastic representation: `U^beta` is Gamma with shape
`1 + 1/beta`, a uniform variate is drawn conditionally, and the product
is wrapped when `wrap = TRUE`.

## Usage

``` r
rsgn(n, mu = 0, sigma = 1, beta = 2, lam1 = 0.5, lam2 = 0.5, wrap = FALSE)
```

## Arguments

- n:

  sample size.

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

- wrap:

  logical; reduce modulo `2*pi`.

## Value

a numeric vector of draws.

## See also

[`rwsgn`](https://mtloots.github.io/wsgnstat/reference/rwsgn.md) for
draws from the wrapped law in its identified parameterisation.

## Examples

``` r
set.seed(2)
x <- rsgn(5000, mu = 0, sigma = 1, beta = 1.3, lam1 = 0.3, lam2 = 0.8)
hist(x, breaks = 60, freq = FALSE, main = "")
curve(dsgn(x, 0, 1, 1.3, 0.3, 0.8), add = TRUE)

## mass left of the location: lam1 / (lam1 + lam2) = 0.2727
mean(x < 0)
#> [1] 0.2832
```
