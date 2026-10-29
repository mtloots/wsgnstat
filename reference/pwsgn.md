# Distribution and quantile functions of the wrapped skew generalised normal distribution

The distribution function is measured anticlockwise from angle `0`:
`pwsgn(theta)` is the probability of the arc `[0, theta]`, computed as
the sum over turns of the line law's mass on `[2 pi k, theta + 2 pi k]`
through the regularised incomplete gamma function, so no numerical
quadrature is involved. The quantile function inverts it by bisection.
The identified parameterisation of
[`dwsgn`](https://mtloots.github.io/wsgnstat/reference/dwsgn.md) is
used.

## Usage

``` r
pwsgn(q, mu = 0, beta = 2, lam = 0.5, K = 50)

qwsgn(p, mu = 0, beta = 2, lam = 0.5, K = 50)
```

## Arguments

- q:

  numeric vector of angles in radians.

- mu:

  location (junction) parameter; for a skew law this is not the mean
  direction.

- beta:

  shape (tail) parameter (\> 0).

- lam:

  skewness parameter in (0, 1); `0.5` is symmetric.

- K:

  truncation of the wrapping sum (`abs(k) <= K`). The default `50`
  captures the whole mass to double precision for `beta >= 0.8`; for
  very heavy tails (`beta` near `0.3`) it does not, and
  `pwsgn(2 * pi, mu, beta, lam, K)` reports the mass that the truncated
  sum captures, so `K` can be raised until that is one.

- p:

  numeric vector of probabilities.

## Value

`pwsgn` gives probabilities in `[0, 1]`: `0` for `q <= 0` and, for
`q >= 2 * pi`, the mass captured by the wrapping sum truncated at `K`
(one to double precision unless the tails are very heavy; see
[`dwsgn`](https://mtloots.github.io/wsgnstat/reference/dwsgn.md)).
`qwsgn` gives angles in `[0, 2 * pi]`, with a warning for any `p` above
the captured mass.

## Examples

``` r
## the distribution function is the integral of the density
pwsgn(2, mu = 1, beta = 1.5, lam = 0.3)
#> [1] 0.6921913
integrate(dwsgn, 0, 2, mu = 1, beta = 1.5, lam = 0.3)$value
#> [1] 0.6921913
## the quantile function inverts it
p <- c(0.05, 0.25, 0.5, 0.75, 0.95)
q <- qwsgn(p, mu = 1, beta = 1.5, lam = 0.3)
rbind(p, pwsgn(q, mu = 1, beta = 1.5, lam = 0.3))
#>   [,1] [,2] [,3] [,4] [,5]
#> p 0.05 0.25  0.5 0.75 0.95
#>   0.05 0.25  0.5 0.75 0.95
## the curve rises from 0 to 1 around the circle
th <- seq(0, 2 * pi, length.out = 361)
plot(th, pwsgn(th, mu = 1, beta = 1.5, lam = 0.3), type = "l",
     xlab = "angle (radians)", ylab = "distribution function")
```
