# wsgnstat 0.1.0

* Initial CRAN submission.
* Density, distribution function, characteristic function and random generation for the skew
  generalised normal distribution on the line (`dsgn()`, `psgn()`, `cfsgn()`, `rsgn()`).
* Density, distribution function, quantile function and random generation for its wrapped version
  on the circle in the identified parameterisation (`dwsgn()`, `pwsgn()`, `qwsgn()`, `rwsgn()`),
  and its trigonometric moments (`trig_moments()`).
* Estimation by minimum distance on the trigonometric moments with sandwich standard errors
  (`fit_wsgn_trig()`) and by maximum likelihood with parametric-bootstrap standard errors
  (`fit_wsgn_mle()`).
* The characteristic function is evaluated through absolutely convergent integrals, so every
  routine is available for all shape parameters `beta > 0`.
* The computation runs on a pure-C back-end shared with the Python package 'wsgnstat'.
