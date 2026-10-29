# Characteristic function of the skew generalised normal distribution (C back-end)

Evaluated through the absolutely convergent cosine and sine integrals
underlying the Fox H-function representation, hence valid for every
`beta > 0`.

## Usage

``` r
cfsgn(t, mu = 0, sigma = 1, beta = 2, lam1 = 0.5, lam2 = 0.5)
```

## Arguments

- t:

  numeric vector of arguments.

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

a complex vector of the same length as `t`.

## Examples

``` r
## reached through absolutely convergent integrals, so every beta > 0 is available,
## including beta <= 1, where the Fox-Wright series does not converge for all t
cfsgn(0, 0, 1, 2, 0.5, 0.5)
#> [1] 1+0i
## the normal member: exp(-t^2 / 4) for N(0, 1/2)
t <- c(0.5, 1.5, 3)
rbind(engine = Re(cfsgn(t, 0, 0.5, 2, 0.5, 0.5)), closed_form = exp(-t^2 / 4))
#>                  [,1]      [,2]      [,3]
#> engine      0.9394131 0.5697828 0.1053992
#> closed_form 0.9394131 0.5697828 0.1053992
## a skew member with beta = 0.8, against the Fourier transform of the density
x <- seq(-60, 60, length.out = 120001)
f <- dsgn(x, 0, 1, 0.8, 0.8, 0.3)
c(engine = cfsgn(1.1, 0, 1, 0.8, 0.8, 0.3),
  quadrature = sum(exp(1i * 1.1 * x) * f) * (x[2] - x[1]))
#>                engine            quadrature 
#> 0.1129021-0.04099506i 0.1129029-0.04100219i
```
