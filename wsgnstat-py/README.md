# wsgnstat (Python)

The two-piece exponential power distribution on the line and its wrapped version on the circle,
treated as a uniform scale mixture with generalised gamma mixing. The name follows the authors'
earlier usage, the wrapped skew generalised normal distribution. The numerics run on a pure-C
back-end shared with the R package 'wsgnstat': the same C source is compiled into both, and the two
fronts return byte-identical values for the same inputs.

Functions:

* line distribution: `dsgn`, `psgn`, `rsgn`, `cfsgn` (characteristic function, computed for every
  shape `beta > 0`; `method="series"` sums the convergent series as a check);
* wrapped distribution with parameters `(mu, sigma, lam, beta)`: `dwsgn`, `pwsgn`, `qwsgn`,
  `rwsgn`, `trig_moments`, `trig_decay`, `nll` (omitting `sigma` gives the reciprocal-rate subfamily
  of version 0.1.0);
* maximum likelihood: `fit_wsgn` (four-parameter fit under a shape floor, or the grouped likelihood
  for rounded data via `resolution`; standard errors from the expected information or a parametric
  bootstrap), `profile_wsgn`, `grouped_loglik`;
* the version 0.1.0 subfamily fitters `fit_wsgn_trig` and `fit_wsgn_mle`.

Pure standard library; no NumPy needed.

## Installation

```
pip install wsgnstat
```

## Example

```python
import math
import random
import wsgnstat as w

rng = random.Random(1)
theta = w.rwsgn(500, mu=1.0, beta=1.6, lam=0.3, sigma=0.4, rng=rng)
fit = w.fit_wsgn(theta)
print(fit["par"], fit["se"], fit["se_method"])

# the characteristic function equals one at zero for every shape
print([abs(w.cfsgn([0.0], 0, 0.4, b, 0.3, 0.7)[0] - 1) for b in (0.1, 0.3, 1, 2)])

# rounded data: the grouped likelihood over 10-degree cells
h = math.pi / 18
theta10 = [(round(v / h) * h) % (2 * math.pi) for v in theta]
print(w.fit_wsgn(theta10, resolution=h, se="none")["par"])
```

## Speed and threads

A fit is a search from many starting values, a screen of candidate locations and, at the smallest
shapes, a spike search, and these trials are independent. `fit_wsgn(..., ncores=4)` and
`profile_wsgn(..., ncores=4)` run them on four threads and give a result that is identical to the
last bit for every `ncores` (the best trial is chosen afterwards in a fixed order); the default is
`ncores=1`. Threads need a library built with OpenMP, which `wsgnstat.HAS_OPENMP` reports.

* Linux wheels are built with OpenMP (libgomp is bundled into the wheel).
* macOS wheels are serial, and `ncores` then has no effect. A wheel that bundles libomp collides
  with the copy of libomp that scikit-learn or PyTorch load into the same process, and delocate
  refuses a libomp built for a newer macOS than the wheel targets, so the robust choice is a
  serial wheel. To get threads on macOS, install Homebrew's `libomp` and build from source with
  `WSGN_OPENMP=1 pip install --no-binary wsgnstat wsgnstat`.
* A source build tries OpenMP on Linux and falls back to a serial library, without failing, when
  the compiler has none; `WSGN_OPENMP=0` forces a serial build.

## Licence

GPL-3.
