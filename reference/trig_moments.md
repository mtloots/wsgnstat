# Population trigonometric moments of the wrapped family (C back-end)

Population trigonometric moments of the wrapped family (C back-end)

## Usage

``` r
trig_moments(P = 2, mu = 0, beta = 2, lam = 0.5)
```

## Arguments

- P:

  highest order p.

- mu:

  location (junction) parameter; for a skew law this is not the mean
  direction.

- beta:

  shape (tail) parameter (\> 0).

- lam:

  skewness parameter in (0, 1); `0.5` is symmetric.

## Value

a list with components `A` and `B`, the cosine and sine moments for
`p = 1, ..., P`.

## Examples

``` r
## the reduction to the wrapped normal in the trigonometric moments: exp(-p^2 / 4)
tm <- trig_moments(4, mu = 0, beta = 2, lam = 0.5)
rbind(engine = tm$A, closed_form = exp(-(1:4)^2 / 4))
#>                  [,1]      [,2]      [,3]       [,4]
#> engine      0.7788008 0.3678794 0.1053992 0.01831564
#> closed_form 0.7788008 0.3678794 0.1053992 0.01831564
## a skew member: the first moment against a large simulated sample
tm1 <- trig_moments(1, mu = 1, beta = 1.5, lam = 0.3)
set.seed(1)
th <- rwsgn(1e5, mu = 1, beta = 1.5, lam = 0.3)
rbind(population = c(tm1$A, tm1$B), sample = c(mean(cos(th)), mean(sin(th))))
#>                  [,1]      [,2]
#> population 0.07582601 0.6411281
#> sample     0.07566699 0.6420514
```
