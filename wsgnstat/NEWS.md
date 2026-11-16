# wsgnstat 0.1.1

## Bug fixes

* `cfsgn()` and `trig_moments()` were inaccurate for shapes below about 0.6: the Simpson rule on
  `[0, 37^(1/beta)]`, capped at two million panels, could not resolve the long range, and
  `cfsgn(0)` was 40624 instead of 1 at `beta = 0.1`. The characteristic function is now evaluated
  on a ray in the upper half-plane by the trapezoid rule in log-radius, which removes both the
  oscillation and the long range, for every `beta > 0`. Against 120-digit values of the convergent
  series at `t = 1, 2, 3`, `lam = 0.3`, two scales (0.4 and `sqrt(0.21)`) and shapes 0.1 to 3, the
  relative error is below `1e-14`. At large `|t|` the two halves of the characteristic function
  cancel, so its error there is absolute, of order `tol`: about `1e-15` to `3e-15` for shapes up to about 3, a few
  times `1e-15` for shapes 5 to 20, and up to about `1e-14` for shapes up to 50, at the default
  `tol`, as the help page now says. The tolerance
  is a new argument `tol` (default `1e-12`), relative to the size of each half, and
  `method = "series"` sums the convergent descending (`beta < 1`) or ascending (`beta > 1`)
  series as a check, refusing (`NaN` with a warning) a half whose estimated rounding error exceeds
  `tol`.

* `dwsgn()`, `pwsgn()` and `qwsgn()` truncated the image sum at `K = 50` images, which was
  inaccurate for small shapes (relative error 0.14 at `beta = 0.3`, and most of the mass missing
  below). The remainder beyond the images summed directly is now added by the Euler-Maclaurin
  formula with the exact incomplete-gamma integral. Against 40-digit image sums at six angles,
  `lam = 0.3`, two scales (0.4 and `sqrt(0.21)`) and shapes 0.1, 0.2, 0.3, 0.5, 1 and 2, the
  relative error of `dwsgn()` is below `1e-14` with the default `K = 50`; `K` is kept as an
  argument, and the result does not depend on it. The tail is kept in log space, so shapes far
  below 0.01 are covered (the density is correct to `1e-13` at `beta = 0.001`). `pwsgn(2 * pi)` is
  now one and `qwsgn()` no longer warns about missing mass. Short arcs keep their relative
  accuracy, and the incomplete gamma functions work from the logarithm of their argument, so small
  probabilities and their quantiles are accurate down to the underflow limit at every shape, also
  next to the location (for instance `qwsgn(1e-300)`).

* `cfsgn()` and `trig_moments()` validate `tol` (one finite number in (0, 1)) and raise values below
  `1e-14`, beyond the rounding of the quadrature, to `1e-14` with a warning. Log-likelihoods at shapes below about 0.6, including those of `fit_wsgn_mle()`, change
  accordingly.

* The help pages called the reciprocal-rate subfamily (`sigma = sqrt(lam (1 - lam))`) "the
  identified parameterisation". That was wrong: the full four-parameter wrapped distribution is
  identifiable, and the subfamily is a three-parameter restriction of it. The documentation now
  says so.

* The wrapped log-density put an observation within half a rounding step (about `4e-16`) of the
  location on the wrong side of it, scoring it as sitting on the location. The signed distance is
  now reduced exactly (`remainder()`), so the side is read from its sign. This changed only
  likelihoods in which an observation lies that close to the fitted location, for instance
  `profile_wsgn(1 + (-3:3) * 1e-14, 2, beta_floor = 0.1)`, which is now 211.601 (was 211.673), and a
  fit whose location is within rounding of an observation: the profile at the fit's own shape now
  equals the fit. Fits on ordinary data are unchanged. The Nelder-Mead search now also stops when
  its simplex has collapsed to a few units in the last place in every coordinate, which the
  corrected likelihood next to an observation would otherwise leave running to its evaluation cap
  (three to four times the run time at shapes below one).

* `fit_wsgn()`: `at_one_sided` now also applies to grouped fits (judged from the log-rates, as for
  exact data), with the same message and the same refusal of `se = "information"`, which now names
  the one-sided case. `at_upper` is flagged within a relative `1e-5` of the upper bound (was
  `1e-6`). `ngrid` is limited to `1e8`, with an error above it. The help page now says that the
  `lam` standard error of a symmetric fit is `NaN` (`lam` is fixed), and that the 1.25 and 1.6
  spike shapes are skipped when they reach half the upper bound (floor 6.25 or more).

* A one-sided fit could stop short of the maximum. When the fitted location sits on an observation
  at the cusp of a one-sided likelihood, the final free polishes collapsed along the location
  before the other parameters had moved: with `set.seed(101); x <- rwsgn(40, 1, 3, 0.5, sigma =
  0.01); th <- signif(1 + abs(x - 1), 15)`, `fit_wsgn(th)` returned 150.39405 at shape 19.68
  while `profile_wsgn(th, 20)` gives 150.40293. The fit now ends with one more polish for fits
  whose log-rates differ by more than `log(1e12)`: the location is held at its fitted value, the
  larger log-rate is pulled back to the smaller plus 45, and the result is kept only if the
  log-likelihood is higher by more than `1e-9`. A free-shape fit that ends on the upper bound of
  the shape is also compared with the search at that bound (the one `profile_wsgn(th, 20)` runs)
  and the higher of the two is kept, so that the result does not depend on whether the compiler
  fuses multiply-adds; `nfev` of such fits rises. The reproducer now gives 150.4029337676 at
  shape 20 (`at_upper` and `at_one_sided` are `TRUE`) on every platform tested. Of 26 test
  sets (20 simulated, the periwinkles at five floors and the ants), compared to eight decimals,
  one other fit changed: a one-sided data set of the same kind rose by 9e-7 (to its profile
  value at shape 20); the other 25 did not.

## New features

* `fit_wsgn()` and `profile_wsgn()` run the independent trials inside a fit in parallel on
  `ncores` threads (`profile_wsgn()` gains the argument; the default is 1 in both, as CRAN
  requires): the starting values, the screen of candidate locations, the spike branches at the
  smallest shapes and their free polishes, and, in `profile_wsgn()`, the shapes. The package must
  be built with OpenMP (4.5 or later); otherwise `ncores` is accepted and ignored. Each trial is
  a function of its own inputs, writes only to its own slot, and the best is chosen afterwards in
  index order (ties to the lowest index), so every result, including the evaluation count, is
  identical to the last bit for every number of threads, and the same as the serial fit. On a
  laptop with four performance and six efficiency cores, busy with other work, a fit of 1500
  observations took 160 s on one thread, 86 s on two, 57 s on four and 55 s on eight (400
  observations: 36, 18, 11 and 13 s); the gain stops near four threads because the first pass of
  the location search, the Nelder-Mead runs and the final polish still run in sequence and the
  extra threads run on slower cores. A likelihood evaluation with 128 or more observations (or
  cells) is split into chunks of 64 over the threads and summed in index order, so it stays
  bit-identical. `ncores` above 256 is reduced to 256 with a message, and the log-gamma function
  is called in its reentrant form (`lgamma_r`) where available, so that concurrent fits share no
  state. Inside
  the parallel bootstrap each refit stays on its own thread. The Python package has the same
  argument and a flag `wsgnstat.HAS_OPENMP`; its Linux wheels are built with OpenMP, its macOS
  wheels are serial (see its README).

* `dwsgn()`, `pwsgn()`, `qwsgn()`, `rwsgn()` and `trig_moments()` take the scale `sigma` as an
  optional argument, giving the full four-parameter wrapped distribution `(mu, sigma, lam, beta)`;
  without it they use the reciprocal-rate subfamily as before. `dwsgn()` gains `log`.

* `fit_wsgn()` fits the four-parameter wrapped distribution by maximum likelihood over
  `beta_floor <= beta <= 20` (default floor 0.1), searching both the regular branch, with the
  location over the observations when the shape is at most one, and the spike branch at three
  shapes, the floor times 1, 1.25 and 1.6, with the location over every observation, and
  returning the best; or by the grouped
  likelihood for rounded data (`resolution`). The regular branch also has one-sided starts, and
  every location search screens each candidate observation after a short optimisation of the
  other parameters (all candidates when there are at most 200, otherwise a structured subset of
  them, as the help page describes); `profile_wsgn()` uses the same search, so the profile at the floor equals the
  fit's floor search. Fits on the shape floor, the upper bound and the rate bound are flagged
  (`at_floor`, `at_upper`, `at_rate_bound`), and the exact log-rates are returned (`log_rates`),
  since `lam` can round to 0 or 1 on a degenerate fit. With `se = "auto"`, an interior fit with
  shape above one gets standard errors from the expected Fisher information (relative difference
  steps, a grid concentrated at the location); otherwise no standard errors are returned, with a
  message saying why. A parametric bootstrap, drawn from the fitted rates, runs only with
  `se = "bootstrap"`. Degenerate data (fewer than three distinct angles, or one occupied cell) and
  non-finite angles are refused with an error. Every step, the optimiser and the bootstrap refits included, runs in C; the refits run
  in parallel when the package is built with OpenMP. `print()`, `coef()`, `vcov()` and `logLik()`
  methods are provided.

* `profile_wsgn()` gives the profile log-likelihood in the shape, `grouped_loglik()` the grouped
  log-likelihood of rounded data, and `trig_decay()` the constant of the exact decay rate
  `C p^-(1 + beta)` of the trigonometric moments.

* `fit_wsgn_trig()` and `fit_wsgn_mle()` are kept for backward compatibility; they fit the
  reciprocal-rate subfamily.

## Documentation

* The example of `fit_wsgn()` now uses 40 observations, and its multi-thread comparison is in
  `\donttest`, to keep its CPU time well under the CRAN limit.

* The title and description now describe the distribution as the wrapped two-piece exponential
  power distribution, treated as a uniform scale mixture with generalised gamma mixing; the package
  name follows the authors' earlier usage, wrapped skew generalised normal.

* The vignette is rewritten around the mixture representation, the characteristic function for
  small and large shapes, the accuracy of the wrapped density, the four-parameter fit with its
  shape floor, and the grouped likelihood for rounded data.

# wsgnstat 0.1.0

* Initial CRAN submission.
* Density, distribution function, characteristic function and random generation for the skew
  generalised normal distribution on the line (`dsgn()`, `psgn()`, `cfsgn()`, `rsgn()`).
* Density, distribution function, quantile function and random generation for its wrapped version
  on the circle in the reciprocal-rate subfamily (`dwsgn()`, `pwsgn()`, `qwsgn()`, `rwsgn()`), and
  its trigonometric moments (`trig_moments()`).
* Estimation by minimum distance on the trigonometric moments with sandwich standard errors
  (`fit_wsgn_trig()`) and by maximum likelihood with parametric-bootstrap standard errors
  (`fit_wsgn_mle()`).
* The computation runs on a pure-C back-end shared with the Python package 'wsgnstat'.
