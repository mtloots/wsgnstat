# wsgnstat

[![R-CMD-check](https://github.com/mtloots/wsgnstat/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mtloots/wsgnstat/actions/workflows/R-CMD-check.yaml)

The **two-piece exponential power distribution** on the line and its **wrapped version on the
circle**, treated as a uniform scale mixture with generalised gamma mixing: given the mixing
variable, the variable is uniform on an interval that contains the location. The package name keeps
the authors' earlier name for the distribution, the wrapped skew generalised normal.

The characteristic function and the trigonometric moments are computed for every shape parameter,
on a rotated path of integration; the wrapped density adds an Euler-Maclaurin tail to the image
sum, so small shapes are covered; and the four-parameter
wrapped distribution is fitted by maximum likelihood under a lower bound on the shape, or by the
grouped likelihood for rounded data.

One shared pure-C back-end serves two front ends:

	wsgnstat/      R package      (.C bindings)
	wsgnstat-py/   Python package (ctypes; pip builds the same C at install time, and a source
	               checkout compiles it on first import)

The two are checked against each other value by value with tools/parity/parity_run.sh, which first
runs parity_sync_guard.sh to confirm that the two copies of the C source are identical. The
reference values used by the tests are made by the scripts in tools/.
