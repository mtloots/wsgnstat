# wsgnstat

[![R-CMD-check](https://github.com/mtloots/wsgnstat/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mtloots/wsgnstat/actions/workflows/R-CMD-check.yaml)

The **skew generalised normal** distribution on the line and its **wrapping onto the circle**.
Directional data — wind direction, animal headings, phase — need circular families, and the
workhorses there are symmetric or only mildly skew. This family carries a shape parameter for tail
weight and separate weights for asymmetry.

The point of the implementation is the **characteristic function**. Estimating by matching
trigonometric moments needs it, and the usual Fox-Wright series converges only over part of the
parameter space; here it is evaluated through absolutely convergent integrals instead, so every
shape parameter is reachable rather than only the convergent region.

A closed form anchors the whole construction: at shape 2 with equal weights the circular law *is*
the wrapped normal with sigma = 1/sqrt(2), and its trigonometric moments are exactly exp(-p^2/4).

One shared pure-C back-end serves two front ends:

	wsgnstat/      R package      (.C bindings; R CMD check: Status OK)
	wsgnstat-py/   Python package (ctypes; compiles the same C on first import)

The two are checked against each other value by value.
