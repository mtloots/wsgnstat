## Submission

This is the first submission of wsgnstat, at version 0.1.0.

## Test environments

* local macOS 27 (arm64), R 4.6.1, `R CMD check --as-cran`: 0 errors, 0 warnings, 1 note.
* win-builder R-devel and R-release: the same tarball was uploaded on 29 October 2026.

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission. Expected for a first submission. In the local check the same note also listed
  the GitHub URLs (repository, issues page and Actions badge) as returning 404, because the
  repository was private when the check was run; it has been made public for this submission and
  the URLs resolve.

## Reverse dependencies

None; this is a new package.

## Notes for the reviewer

The package compiles a small pure-C back-end (no OpenMP, no random numbers in C, no printing),
registered through `R_registerRoutines` and called with `.C`. The same C source also backs a
Python package of the same name, and the two return byte-identical values for the same inputs.

The correctness claim is that the densities, distribution functions, quantiles, characteristic
function and trigonometric moments are right for every shape parameter, and that the estimators
recover the parameters that generated the data. The tests establish this against known values, not
merely by running the code:

* the wrapped density equals the wrapped normal density (sigma = 1/sqrt(2)) at beta = 2 and
  lam = 0.5 to 1e-9, and the trigonometric moments equal exp(-p^2/4) to 1e-9;
* the line and wrapped densities integrate to one, and the line and wrapped distribution functions
  agree with numerical integration of the densities to 1e-9, with the normal distribution function
  and the wrapped normal distribution function as closed-form anchors;
* the characteristic function agrees with the Fourier transform of the density computed by
  quadrature, including at a skew, heavy-tailed member;
* the quantile function inverts the distribution function to 1e-12;
* draws from `rwsgn()` pass a Kolmogorov-Smirnov test against the exact distribution function at
  the 0.1 per cent level;
* for three known parameter sets (including beta = 0.8, below the region where the usual series for
  the characteristic function applies) the moment estimator lands within four of its sandwich
  standard errors of the truth and the maximum likelihood estimator close to it.

Examples run in under two seconds each; one bootstrap example is wrapped in `\donttest{}` because it
refits the model once per replicate. The vignette uses the 'circular' package (Suggests) for a
wind-direction data set and is skipped gracefully when it is unavailable.
