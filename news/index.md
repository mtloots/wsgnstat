# Changelog

## wsgnstat 0.1.0

- Initial CRAN submission.
- Density, distribution function, characteristic function and random
  generation for the skew generalised normal distribution on the line
  ([`dsgn()`](https://mtloots.github.io/wsgnstat/reference/dsgn.md),
  [`psgn()`](https://mtloots.github.io/wsgnstat/reference/psgn.md),
  [`cfsgn()`](https://mtloots.github.io/wsgnstat/reference/cfsgn.md),
  [`rsgn()`](https://mtloots.github.io/wsgnstat/reference/rsgn.md)).
- Density, distribution function, quantile function and random
  generation for its wrapped version on the circle in the identified
  parameterisation
  ([`dwsgn()`](https://mtloots.github.io/wsgnstat/reference/dwsgn.md),
  [`pwsgn()`](https://mtloots.github.io/wsgnstat/reference/pwsgn.md),
  [`qwsgn()`](https://mtloots.github.io/wsgnstat/reference/pwsgn.md),
  [`rwsgn()`](https://mtloots.github.io/wsgnstat/reference/rwsgn.md)),
  and its trigonometric moments
  ([`trig_moments()`](https://mtloots.github.io/wsgnstat/reference/trig_moments.md)).
- Estimation by minimum distance on the trigonometric moments with
  sandwich standard errors
  ([`fit_wsgn_trig()`](https://mtloots.github.io/wsgnstat/reference/fit_wsgn_trig.md))
  and by maximum likelihood with parametric-bootstrap standard errors
  ([`fit_wsgn_mle()`](https://mtloots.github.io/wsgnstat/reference/fit_wsgn_mle.md)).
- The characteristic function is evaluated through absolutely convergent
  integrals, so every routine is available for all shape parameters
  `beta > 0`.
- The computation runs on a pure-C back-end shared with the Python
  package ‘wsgnstat’.
