# wsgnstat (Python)

The skew generalised normal distribution on the line and its wrapped version on the circle, on a
pure-C back-end shared with the R package 'wsgnstat': the same C source is compiled into both, and
the two fronts return the same values for the same inputs.

Functions: `dsgn`, `psgn`, `cfsgn`, `rsgn` (line law); `dwsgn`, `pwsgn`, `qwsgn`, `rwsgn` (wrapped
law in the identified parameterisation: location `mu`, tail `beta`, skewness `lam`);
`trig_moments`; `nll`; `fit_wsgn_trig` (minimum distance on trigonometric moments, sandwich
standard errors) and `fit_wsgn_mle` (maximum likelihood, optional parametric-bootstrap standard
errors). The characteristic function is evaluated through absolutely convergent integrals, so every
shape parameter `beta > 0` is available. Pure standard library; no NumPy needed.

## Installation

```
pip install wsgnstat
```

## Example

```python
import random
import wsgnstat as w

rng = random.Random(1)
theta = w.rwsgn(1000, mu=1.0, beta=1.5, lam=0.3, rng=rng)
fit = w.fit_wsgn_trig(theta)
print(fit["par"], fit["se"])

# at beta = 2 and lam = 0.5 the circular law is the wrapped normal with sigma = 1/sqrt(2),
# whose trigonometric moments are exp(-p^2/4)
print(w.trig_moments(3, 0.0, 2.0, 0.5)["A"])
```

## Licence

GPL-3.
