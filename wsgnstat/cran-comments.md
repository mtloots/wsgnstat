## Submission

This is an update of wsgnstat from 0.1.0 (on CRAN since 14 November 2026) to 0.1.1. It fixes
numerical defects found after the release and adds the four-parameter likelihood:

* `cfsgn()` and `trig_moments()` were inaccurate for shapes below about 0.6 (the characteristic
  function at zero was 40624 instead of 1 at shape 0.1), because Simpson's rule on a very long
  range could not resolve it. The integral is now taken along a rotated path by an exponentially
  convergent trapezoid rule for every shape; against 120-digit series values for shapes 0.1 to 3
  (t = 1, 2, 3, two scales) the relative error is below 1e-14.
* `dwsgn()`, `pwsgn()` and `qwsgn()` truncated the image sum at 50 images, which was inaccurate for
  small shapes (relative error 0.14 at shape 0.3). An Euler-Maclaurin tail with the exact
  incomplete-gamma integral now corrects this; against 40-digit image sums for shapes 0.1 to 2
  (six angles, two scales) the relative error of the density is below 1e-14.
* `fit_wsgn()` and `profile_wsgn()` take `ncores` (default 1) and run the independent trials inside
  a fit on that many OpenMP threads, with a result that is bit-identical for every number of
  threads; the package builds and runs without OpenMP (via `$(SHLIB_OPENMP_CFLAGS)`), where
  `ncores` is ignored. Examples and tests use at most two threads.
* A one-sided fit that stopped short of the maximum at a cusp is polished further.
* The help pages wrongly called a three-parameter subfamily "the identified parameterisation"; the
  full four-parameter family is identifiable, and the documentation now says so.
* New: the scale argument `sigma` for the wrapped functions, `fit_wsgn()` (maximum likelihood for
  the four-parameter wrapped distribution, with a shape floor, the grouped likelihood for rounded
  data, and information or bootstrap standard errors), `profile_wsgn()`, `grouped_loglik()` and
  `trig_decay()`. The fit covers both the regular branch and the spike branch at the shape floor,
  flags fits on any bound, and refuses degenerate or non-finite data.
* The Title and Description now describe the distribution under its usual name, the wrapped
  two-piece exponential power distribution; the package name is unchanged.

## Test environments

* local macOS 27 (arm64), R 4.6.1, `R CMD check --as-cran` on the built tarball, with
  `_R_CHECK_FORCE_SUGGESTS_=true`.
* Debian (x86_64 and aarch64), R 4.6.2, gcc 16 with OpenMP: the full test suite, and on aarch64
  `R CMD check --as-cran`.
* Windows (MinGW-w64, msvcrt and UCRT): the C sources cross-compiled and linked, with and without
  OpenMP.
* GitHub Actions R-CMD-check matrix (Ubuntu R devel, release and oldrel-1; macOS and Windows R release).

## R CMD check results

0 errors | 0 warnings | 1 note

* "Days since last update: 2". Version 0.1.0 was published on 14 November 2026; this update follows
  quickly because it corrects numerical errors in the characteristic function and the wrapped
  density at small shape parameters, which users of 0.1.0 would otherwise meet.

## Reverse dependencies

There are no reverse dependencies.

## Notes for the reviewer

The package compiles a pure-C back-end registered through `R_registerRoutines` and called with
`.C`; it uses OpenMP, through `SHLIB_OPENMP_CFLAGS`, for the independent trials inside a fit or
profile and for the bootstrap refits, on `ncores` threads (default 1), with the same result for
every number of threads. The C code draws no random
numbers from R: bootstrap samples are drawn in R with `set.seed()` reproducibility, and the only
random stream in C is a fixed-seed generator for starting values. The same C source backs a Python
package of the same name, and the two return byte-identical values.

The correctness tests compare against known values: the characteristic function against
120-digit values of its convergent series for shapes 0.1 to 3 and against closed forms at shapes
one and two, the wrapped density against 40-digit Euler-Maclaurin image sums, the distribution
functions against numerical integration, the mixture representation by simulation, and the
maximum likelihood fits against values computed independently for two data sets of the 'circular'
package. Examples run in a few seconds; the bootstrap example is wrapped in `\donttest{}`.
