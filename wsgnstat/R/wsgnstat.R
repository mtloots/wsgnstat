## wsgnstat: the two-piece exponential power distribution on the line and its wrapped version on
## the circle, treated as a uniform scale mixture with generalised gamma mixing. The package keeps
## the authors' earlier name for this distribution, the (wrapped) skew generalised normal. All the
## computation, including every loop, optimiser and bootstrap refit, runs in one shared pure-C
## back-end compiled at install, which also backs the Python package 'wsgnstat'.
#' @useDynLib wsgnstat, .registration = TRUE, .fixes = "C_"
#' @importFrom stats optim rgamma runif sd cov setNames coef vcov logLik
#' @keywords internal
"_PACKAGE"

## ---------------------------------------------------------------- internal helpers

## the scale of the reciprocal-rate subfamily (ab = 1) when sigma is not given
.sigma_or_default <- function(sigma, lam) {
  if (is.null(sigma)) sqrt(pmax(lam * (1 - lam), 0)) else sigma
}

.check_line <- function(sigma, beta, lam1, lam2) {
  ok <- function(v) length(v) == 1L && is.numeric(v) && is.finite(v) && v > 0
  if (!ok(sigma) || !ok(beta) || !ok(lam1) || !ok(lam2))
    stop("'sigma', 'beta', 'lam1' and 'lam2' must each be one positive finite number")
  invisible(TRUE)
}

.check_tol <- function(tol) {
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0 || tol >= 1)
    stop("'tol' must be one finite number in (0, 1)")
  if (tol < 1e-14) {
    warning("'tol' below 1e-14 is beyond the rounding of the quadrature; 1e-14 is used")
    tol <- 1e-14
  }
  tol
}

.check_K <- function(K) {
  if (length(K) != 1L || !is.numeric(K) || !is.finite(K) || K < 1 || K > 1e6 || K != round(K))
    stop("'K' must be a whole number between 1 and 1e6")
  invisible(TRUE)
}

.check_wrapped <- function(sigma, beta, lam, K = 1L) {
  if (length(lam) != 1L || !is.finite(lam) || lam <= 0 || lam >= 1) stop("'lam' must lie in (0, 1)")
  if (length(sigma) != 1L || !is.finite(sigma) || sigma <= 0)
    stop("'sigma' must be one positive number")
  if (length(beta) != 1L || !is.finite(beta) || beta <= 0)
    stop("'beta' must be one positive number")
  .check_K(K)
  invisible(TRUE)
}

#' Density of the two-piece exponential power distribution on the line (C back-end)
#'
#' The density is \code{c exp(-(lam1 (x - mu) / sigma)^beta)} for \code{x >= mu} and
#' \code{c exp(-(lam2 (mu - x) / sigma)^beta)} for \code{x < mu}. With \code{lam1 = lam} and
#' \code{lam2 = 1 - lam} it is the distribution of \code{mu + Z V}, where \code{Z^beta} is gamma
#' distributed with shape \code{1 + 1/beta} and \code{V} is uniform on
#' \code{[-sigma / (1 - lam), sigma / lam]}, independent of \code{Z}: a uniform scale mixture with
#' generalised gamma mixing, with \code{P(X < mu) = lam}.
#' @param x numeric vector of quantiles.
#' @param mu location parameter.
#' @param sigma scale parameter (> 0).
#' @param beta shape parameter (> 0).
#' @param lam1 right-hand rate parameter (> 0); the right-hand rate is \code{lam1 / sigma}.
#' @param lam2 left-hand rate parameter (> 0); the left-hand rate is \code{lam2 / sigma}.
#' @return the probability density function evaluated at \code{x}.
#' @examples
#' ## a right-skew member with a cusp at its location, and the normal member
#' x <- seq(-4, 6, length.out = 401)
#' f_skew <- dsgn(x, mu = 0, sigma = 1, beta = 0.8, lam1 = 0.3, lam2 = 0.7)
#' f_norm <- dsgn(x, mu = 0, sigma = 0.5, beta = 2, lam1 = 0.5, lam2 = 0.5)
#' plot(x, f_skew, type = "l", ylab = "density")
#' lines(x, f_norm, lty = 2)
#' ## at beta = 2 with lam1 = lam2 = 0.5 and sigma = 0.5 the distribution is N(0, 1/2)
#' max(abs(f_norm - dnorm(x, 0, sqrt(1 / 2))))
#' ## the density integrates to one
#' integrate(dsgn, -Inf, Inf, mu = 0, sigma = 1, beta = 0.8, lam1 = 0.3, lam2 = 0.7)$value
#' @export
dsgn <- function(x, mu = 0, sigma = 1, beta = 2, lam1 = 0.5, lam2 = 0.5) {
  .check_line(sigma, beta, lam1, lam2)
  .C(C_wsgn_density, x = as.double(x), nx = length(x), mu = as.double(mu),
     sigma = as.double(sigma), beta = as.double(beta), lam1 = as.double(lam1),
     lam2 = as.double(lam2), out = double(length(x)), NAOK = TRUE)$out
}

#' Distribution function of the two-piece exponential power distribution on the line
#'
#' Evaluated in closed form through the regularised incomplete gamma function: the mass to the left
#' of \code{mu} is \code{lam1 / (lam1 + lam2)}, and each tail is a scaled gamma tail in
#' \code{(lam * |x - mu| / sigma)^beta}.
#' @param q numeric vector of quantiles.
#' @inheritParams dsgn
#' @return the distribution function evaluated at \code{q}.
#' @examples
#' ## agreement with numerical integration of the density
#' q <- c(-1, 0, 0.7, 2.5)
#' closed <- psgn(q, 0, 1, 1.3, 0.3, 0.8)
#' numeric <- vapply(q, function(u)
#'   integrate(dsgn, -Inf, u, mu = 0, sigma = 1, beta = 1.3, lam1 = 0.3, lam2 = 0.8)$value, 0)
#' rbind(closed, numeric)
#' ## the normal member: psgn equals pnorm with variance 1/2
#' psgn(0.4, 0, 0.5, 2, 0.5, 0.5) - pnorm(0.4, 0, sqrt(1 / 2))
#' @export
psgn <- function(q, mu = 0, sigma = 1, beta = 2, lam1 = 0.5, lam2 = 0.5) {
  .check_line(sigma, beta, lam1, lam2)
  .C(C_wsgn_pdist, x = as.double(q), nx = length(q), mu = as.double(mu),
     sigma = as.double(sigma), beta = as.double(beta), lam1 = as.double(lam1),
     lam2 = as.double(lam2), out = double(length(q)), NAOK = TRUE)$out
}

#' Density of the wrapped two-piece exponential power distribution (C back-end)
#'
#' The wrapped distribution is that of \code{X mod 2 pi}, with \code{X} the line distribution of
#' \code{\link{dsgn}} with right rate \code{a = lam / sigma} and left rate
#' \code{b = (1 - lam) / sigma}. Its four parameters are the location \code{mu}, the scale
#' \code{sigma}, the skewness \code{lam} (\code{P(X < mu) = lam}) and the shape \code{beta}, and
#' the four-parameter family is identifiable. When \code{sigma} is not given, the reciprocal-rate
#' subfamily \code{a b = 1}, that is \code{sigma = sqrt(lam (1 - lam))}, is used, as in version
#' 0.1.0; it is a three-parameter restriction of the model, kept for backward compatibility.
#'
#' The density is the image sum \code{c sum_k [exp(-(a (t + 2 pi k))^beta) + exp(-(b (2 pi - t + 2
#' pi k))^beta)]}, \code{t = (theta - mu) mod 2 pi}, with \code{c = lam (1 - lam) / (sigma Gamma(1 +
#' 1/beta))}. On each side of \code{mu} images are summed directly, and the remainder is added by
#' the Euler-Maclaurin formula, with the exact incomplete-gamma integral and Bernoulli terms through
#' \code{B_12}, as soon as it is accurate; it is tried after 8, 16, 32 and \code{K} images and after
#' each further one, so the result does not depend on \code{K}. The tail is kept in log space, so
#' very small shapes are covered (at \code{beta = 0.005} and \code{0.001} the density differs from
#' \code{1/(2 pi)} by less than \code{1e-13} relative, the limit set by cancelling
#' \code{lgamma(1/beta)}). Against 40-digit image sums at six angles, \code{lam = 0.3}, \code{sigma
#' = 0.4} and \code{sqrt(0.21)}, and shapes 0.1, 0.2, 0.3, 0.5, 1 and 2 (the file
#' \code{extdata/ref_dw.csv}), the relative error is below \code{1e-14}.
#' @param theta numeric vector of angles in radians.
#' @param mu location parameter.
#' @param beta shape parameter (> 0).
#' @param lam skewness parameter in (0, 1); \code{0.5} is symmetric.
#' @param K number of images summed directly on each side of \code{mu} before the
#'   Euler-Maclaurin tail (a positive integer); more are added directly when the formula is not
#'   yet accurate there, and fewer are used when the rest is negligible.
#' @param sigma scale parameter (> 0); \code{NULL} gives the reciprocal-rate subfamily.
#' @param log logical; return the log-density.
#' @return the wrapped density (or its logarithm) evaluated at \code{theta}.
#' @examples
#' ## a skew member with a cusp, drawn on the circle as a radial curve about the unit circle
#' th <- seq(0, 2 * pi, length.out = 721)
#' d <- dwsgn(th, mu = pi / 2, beta = 0.8, lam = 0.3, sigma = 0.4)
#' plot(cos(th), sin(th), type = "l", asp = 1, xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5),
#'      xlab = "", ylab = "", lty = 3)
#' lines((1 + d) * cos(th), (1 + d) * sin(th))
#' ## the density integrates to one over the circle, also at a very small shape
#' integrate(dwsgn, 0, 2 * pi, mu = 1, beta = 0.2, lam = 0.3, sigma = 0.4,
#'           subdivisions = 1000)$value
#' ## at beta = 2 with lam = 0.5 the circular distribution is the wrapped normal
#' ## with standard deviation sqrt(2) sigma
#' tt <- c(0.5, 2.0, 4.5)
#' wn <- vapply(tt, function(t) sum(dnorm(t + 2 * pi * (-80:80), 0, sqrt(2) * 0.6)), 0)
#' rbind(wsgnstat = dwsgn(tt, 0, 2, 0.5, sigma = 0.6), wrapped_normal = wn)
#' @export
dwsgn <- function(theta, mu = 0, beta = 2, lam = 0.5, K = 50, sigma = NULL, log = FALSE) {
  sigma <- .sigma_or_default(sigma, lam)
  .check_wrapped(sigma, beta, lam, K)
  .C(C_wsgn_wdensity, theta = as.double(theta), nt = length(theta), mu = as.double(mu),
     sigma = as.double(sigma), beta = as.double(beta), lam = as.double(lam),
     K = as.integer(K), give_log = as.integer(isTRUE(log)), out = double(length(theta)),
     NAOK = TRUE)$out
}

#' Distribution and quantile functions of the wrapped two-piece exponential power distribution
#'
#' The distribution function is measured anticlockwise from angle \code{0}: \code{pwsgn(theta)} is
#' the probability of the arc \code{[0, theta]}. It is computed from the incomplete-gamma survival
#' functions of the two sides, summed over images, with the Euler-Maclaurin tail beyond the
#' \code{K}-th image. The arc is described by its start and its width, so a short arc keeps its
#' relative accuracy: for an arc shorter than 0.01 the first image is taken from the incomplete
#' gamma function where that is accurate and by Gauss-Legendre quadrature of the density otherwise,
#' and the incomplete gamma functions are evaluated from the logarithm of their argument, so that
#' \code{(rate x)^beta} never underflows next to the location, at any shape. Small probabilities
#' therefore keep their relative accuracy down to the underflow limit, on both sides of the location
#' and for large shapes too (checked for shapes 0.3 to 50). The quantile function inverts it by
#' bisection with geometric midpoints while the bracket is wide, so that small probabilities get
#' quantiles with relative accuracy (for instance \code{qwsgn(1e-300)}). The parameters are as in
#' \code{\link{dwsgn}}.
#' @param q numeric vector of angles in radians.
#' @param p numeric vector of probabilities.
#' @inheritParams dwsgn
#' @return \code{pwsgn} gives probabilities in \code{[0, 1]}: \code{0} for \code{q <= 0} and
#'   \code{1} for \code{q >= 2 * pi}. \code{qwsgn} gives angles in \code{[0, 2 * pi]}.
#' @examples
#' ## the distribution function is the integral of the density
#' pwsgn(2, mu = 1, beta = 1.5, lam = 0.3, sigma = 0.5)
#' integrate(dwsgn, 0, 2, mu = 1, beta = 1.5, lam = 0.3, sigma = 0.5)$value
#' ## the quantile function inverts it
#' p <- c(0.05, 0.25, 0.5, 0.75, 0.95)
#' q <- qwsgn(p, mu = 1, beta = 1.5, lam = 0.3, sigma = 0.5)
#' rbind(p, pwsgn(q, mu = 1, beta = 1.5, lam = 0.3, sigma = 0.5))
#' ## the curve rises from 0 to 1 around the circle
#' th <- seq(0, 2 * pi, length.out = 361)
#' plot(th, pwsgn(th, mu = 1, beta = 1.5, lam = 0.3, sigma = 0.5), type = "l",
#'      xlab = "angle (radians)", ylab = "distribution function")
#' @name pwsgn
NULL

#' @rdname pwsgn
#' @export
pwsgn <- function(q, mu = 0, beta = 2, lam = 0.5, K = 50, sigma = NULL) {
  sigma <- .sigma_or_default(sigma, lam)
  .check_wrapped(sigma, beta, lam, K)
  .C(C_wsgn_wcdf, theta = as.double(q), nt = length(q), mu = as.double(mu),
     sigma = as.double(sigma), beta = as.double(beta), lam = as.double(lam),
     K = as.integer(K), out = double(length(q)), NAOK = TRUE)$out
}

#' @rdname pwsgn
#' @export
qwsgn <- function(p, mu = 0, beta = 2, lam = 0.5, K = 50, sigma = NULL) {
  sigma <- .sigma_or_default(sigma, lam)
  .check_wrapped(sigma, beta, lam, K)
  .C(C_wsgn_wquant, p = as.double(p), np = length(p), mu = as.double(mu),
     sigma = as.double(sigma), beta = as.double(beta), lam = as.double(lam),
     K = as.integer(K), out = double(length(p)), NAOK = TRUE)$out
}

#' Characteristic function of the two-piece exponential power distribution (C back-end)
#'
#' The characteristic function is
#' \code{exp(i t mu) [w1 G(t sigma / lam1) + w2 Conj(G(t sigma / lam2))] / Gamma(1 + 1/beta)},
#' with \code{w1 = lam2 / (lam1 + lam2) = P(X > mu)}, \code{w2 = 1 - w1} and
#' \code{G(w) = int_0^Inf exp(i w u - u^beta) du}. The default method evaluates \code{G} on the
#' ray \code{u = r exp(i phi)} in the upper half-plane, with \code{phi} half of
#' \code{min(pi, pi / (2 beta))}, by the trapezoid rule in \code{log r}, halving the step until
#' successive sums agree to \code{tol}; the rotation removes both the oscillation and the very long
#' range that the integral has on the real axis for small shapes.
#'
#' \emph{Accuracy.} \code{tol} is relative to the size of each half, \code{G(w) / Gamma(1 +
#' 1/beta)}, which is of order \code{min(1, 1/|w|)}. The two halves cancel to leading order at large
#' \code{|t|}, so the error of the characteristic function is absolute, of order \code{tol min(1,
#' max(lam1, lam2) / (sigma |t|))} (in practice of order \code{tol}: about \code{1e-15} to \code{3e-15} for shapes up to about 3, a few times
#' \code{1e-15} for shapes 5 to 20, and up to about \code{1e-14} for shapes up to 50, at the default
#' \code{tol}); values much smaller
#' than that carry no relative accuracy (for the normal member, for instance, beyond \code{t} of
#' about 8 in standard-deviation units). Against 120-digit values of the convergent series at
#' \code{t = 1, 2, 3}, \code{lam1 = 0.3}, \code{lam2 = 0.7}, \code{sigma = 0.4} and
#' \code{sqrt(0.21)}, and shapes 0.1 to 3 (the file \code{extdata/ref_cf.csv}), the relative error
#' is below \code{1e-14}.
#'
#' \emph{Series.} The method \code{"series"} sums the convergent series instead, as a check: the
#' descending series for \code{beta < 1}, the ascending (Fox-Wright) series for \code{beta > 1},
#' and \code{1 / (1 - i w)} at \code{beta = 1}. Its rounding error is estimated as
#' \code{1e-15 sum_k |t_k| (1 + |log t_k| + |arg t_k|)} over the terms \code{t_k}, and a half is
#' refused, with \code{NaN} and a warning, when that estimate exceeds \code{tol} relative to the
#' half, when the series does not converge within 200000 terms, or when the sum breaks the bound
#' \code{|G(w)| <= Gamma(1 + 1/beta)}. On scans of 400 arguments per shape for shapes 0.05 to 50 the
#' accepted values agreed with the quadrature to about \code{tol} relative to each half (to within a
#' factor 3 at \code{tol = 1e-8}); like the quadrature, the combined value is accurate only in
#' absolute terms at large \code{|t|}. Many arguments are refused, more for shapes near one.
#' @param t numeric vector of arguments.
#' @inheritParams dsgn
#' @param tol tolerance relative to the size of each half (see Accuracy), for the quadrature and for
#'   the estimated rounding error of a series: one finite number in (0, 1). Values below
#'   \code{1e-14}, beyond the rounding of the quadrature, are raised to \code{1e-14} with a warning.
#' @param method \code{"quadrature"} (default) or \code{"series"}.
#' @return a complex vector of the same length as \code{t}.
#' @examples
#' ## every shape is reachable: the value at zero is one, also for small shapes
#' sapply(c(0.1, 0.2, 0.3, 0.5, 1, 2), function(b) Mod(cfsgn(0, 0, 0.4, b, 0.3, 0.7) - 1))
#' ## below shape one the descending series converges and agrees with the quadrature
#' t <- c(1, 2, 3)
#' cbind(quadrature = cfsgn(t, 0, 0.4, 0.3, 0.3, 0.7),
#'       series = cfsgn(t, 0, 0.4, 0.3, 0.3, 0.7, method = "series"))
#' ## the normal member: exp(-t^2 / 4) for N(0, 1/2)
#' rbind(engine = Re(cfsgn(t, 0, 0.5, 2, 0.5, 0.5)), closed_form = exp(-t^2 / 4))
#' ## the asymmetric Laplace member (beta = 1) in closed form
#' a <- 0.3 / 0.4
#' b <- 0.7 / 0.4
#' rbind(engine = cfsgn(t, 0, 0.4, 1, 0.3, 0.7),
#'       closed_form = a * b / ((a - 1i * t) * (b + 1i * t)))
#' @export
cfsgn <- function(t, mu = 0, sigma = 1, beta = 2, lam1 = 0.5, lam2 = 0.5, tol = 1e-12,
                  method = c("quadrature", "series")) {
  method <- match.arg(method)
  .check_line(sigma, beta, lam1, lam2)
  tol <- .check_tol(tol)
  r <- .C(C_wsgn_cf, t = as.double(t), nt = length(t), mu = as.double(mu),
          sigma = as.double(sigma), beta = as.double(beta), lam1 = as.double(lam1),
          lam2 = as.double(lam2), tol = as.double(tol), method = as.integer(method == "series"),
          re = double(length(t)), im = double(length(t)), NAOK = TRUE)
  if (method == "series" && anyNA(r$re))
    warning("the series cancels beyond the requested precision at some arguments (NaN returned); ",
            "use the default quadrature there")
  complex(real = r$re, imaginary = r$im)
}

#' Trigonometric moments of the wrapped distribution (C back-end)
#'
#' The trigonometric moments of the wrapped distribution are the values of the characteristic
#' function of the line distribution at the integers, \code{phi_p = E exp(i p Theta) = A_p + i B_p},
#' computed as in \code{\link{cfsgn}}, so their error is absolute, of order \code{tol / p} (in
#' practice of order \code{tol}: about \code{1e-15} to \code{3e-15} for shapes up to about 3, a few times
#' \code{1e-15} for shapes 5 to 20, and up to about \code{1e-14} for shapes up to 50, at the default
#' \code{tol}); moments much smaller than that (the far moments of a smooth
#' member) carry no relative accuracy. \code{trig_decay} gives the constant \code{C} of the exact
#' decay rate \code{|phi_p| ~ C p^-(1 + beta)} as \code{p} grows, which holds unless
#' \code{lam = 0.5} and \code{beta} is an even integer (then \code{C = 0} and the decay is faster
#' than any power).
#' @param P highest order \code{p}, a whole number from 1 to \code{1e7}.
#' @inheritParams dwsgn
#' @param tol tolerance, as in \code{\link{cfsgn}}: one finite number in (0, 1), raised to
#'   \code{1e-14} if smaller.
#' @return \code{trig_moments}: a list with components \code{A} and \code{B}, the cosine and sine
#'   moments for \code{p = 1, ..., P}. \code{trig_decay}: a named vector with the constant \code{C}
#'   and the exponent \code{1 + beta}.
#' @examples
#' ## the reduction to the wrapped normal: exp(-sigma^2 p^2) at beta = 2, lam = 0.5
#' tm <- trig_moments(4, mu = 0, beta = 2, lam = 0.5, sigma = 0.5)
#' rbind(engine = tm$A, closed_form = exp(-0.25 * (1:4)^2))
#' ## a skew member with a small shape: the moments approach C p^-(1 + beta)
#' tm <- trig_moments(400, mu = 0, beta = 0.3, lam = 0.3, sigma = 0.4)
#' dc <- trig_decay(beta = 0.3, lam = 0.3, sigma = 0.4)
#' p <- c(10, 100, 400)
#' Mod(complex(real = tm$A[p], imaginary = tm$B[p])) * p^dc["exponent"] / dc["C"]
#' ## the first moment against a large simulated sample
#' tm1 <- trig_moments(1, mu = 1, beta = 1.5, lam = 0.3, sigma = 0.5)
#' set.seed(1)
#' th <- rwsgn(1e5, mu = 1, beta = 1.5, lam = 0.3, sigma = 0.5)
#' rbind(population = c(tm1$A, tm1$B), sample = c(mean(cos(th)), mean(sin(th))))
#' @export
trig_moments <- function(P = 2, mu = 0, beta = 2, lam = 0.5, sigma = NULL, tol = 1e-12) {
  if (length(P) != 1L || !is.numeric(P) || !is.finite(P) || P < 1 || P != round(P) || P > 1e7)
    stop("'P' must be one whole number between 1 and 1e7")
  sigma <- .sigma_or_default(sigma, lam)
  .check_wrapped(sigma, beta, lam)
  tol <- .check_tol(tol)
  r <- .C(C_wsgn_trigmom, P = as.integer(P), mu = as.double(mu), sigma = as.double(sigma),
          beta = as.double(beta), lam = as.double(lam), tol = as.double(tol),
          A = double(P), B = double(P))
  list(A = r$A, B = r$B)
}

#' @rdname trig_moments
#' @export
trig_decay <- function(beta = 2, lam = 0.5, sigma = NULL) {
  sigma <- .sigma_or_default(sigma, lam)
  .check_wrapped(sigma, beta, lam)
  r <- .C(C_wsgn_decay, sigma = as.double(sigma), beta = as.double(beta), lam = as.double(lam),
          out = double(2))$out
  setNames(r, c("C", "exponent"))
}

#' Random generation from the two-piece exponential power distribution on the line
#'
#' Uses the uniform scale mixture: \code{Z^beta} is drawn from the gamma distribution with shape
#' \code{1 + 1/beta}, and the draw is uniform on \code{[mu - Z sigma / lam2, mu + Z sigma / lam1]},
#' an interval that contains \code{mu}; a single uniform variable decides the side. The draws are
#' reduced modulo \code{2 pi} when \code{wrap = TRUE}.
#' @param n sample size.
#' @inheritParams dsgn
#' @param wrap logical; reduce modulo \code{2 pi}.
#' @return a numeric vector of draws.
#' @seealso \code{\link{rwsgn}} for draws from the wrapped distribution.
#' @examples
#' set.seed(2)
#' x <- rsgn(5000, mu = 0, sigma = 1, beta = 1.3, lam1 = 0.3, lam2 = 0.7)
#' hist(x, breaks = 60, freq = FALSE, main = "")
#' curve(dsgn(x, 0, 1, 1.3, 0.3, 0.7), add = TRUE)
#' ## the mass to the left of the location is lam1 / (lam1 + lam2) = 0.3
#' mean(x < 0)
#' @export
rsgn <- function(n, mu = 0, sigma = 1, beta = 2, lam1 = 0.5, lam2 = 0.5, wrap = FALSE) {
  .check_line(sigma, beta, lam1, lam2)
  u <- rgamma(n, shape = 1 + 1 / beta)^(1 / beta)
  x <- runif(n, min = mu - (sigma / lam2) * u, max = mu + (sigma / lam1) * u)
  if (wrap) x %% (2 * pi) else x
}

#' Random generation from the wrapped two-piece exponential power distribution
#'
#' Draws from the wrapped distribution with the parameters of \code{\link{dwsgn}}, that is
#' \code{rsgn(n, mu, sigma, beta, lam1 = lam, lam2 = 1 - lam, wrap = TRUE)}. When \code{sigma} is
#' not given, the reciprocal-rate subfamily \code{sigma = sqrt(lam (1 - lam))} is used.
#' @inheritParams dwsgn
#' @param n number of draws.
#' @return a numeric vector of angles in \code{[0, 2 * pi)}.
#' @examples
#' set.seed(3)
#' th <- rwsgn(2000, mu = 1, beta = 0.9, lam = 0.3, sigma = 0.4)
#' hist(th, breaks = seq(0, 2 * pi, length.out = 49), freq = FALSE, main = "",
#'      xlab = "angle (radians)")
#' curve(dwsgn(x, mu = 1, beta = 0.9, lam = 0.3, sigma = 0.4), 0, 2 * pi, add = TRUE, n = 721)
#' @export
rwsgn <- function(n, mu = 0, beta = 2, lam = 0.5, sigma = NULL) {
  sigma <- .sigma_or_default(sigma, lam)
  .check_wrapped(sigma, beta, lam)
  rsgn(n, mu = mu, sigma = sigma, beta = beta, lam1 = lam, lam2 = 1 - lam, wrap = TRUE)
}

#' Fit the reciprocal-rate subfamily by trigonometric moments
#'
#' Minimum-distance matching of the first two sample trigonometric moments to those of the
#' three-parameter reciprocal-rate subfamily (\code{sigma = sqrt(lam (1 - lam))}), with standard
#' errors from the minimum-distance sandwich. Kept from version 0.1.0 for backward compatibility.
#' Matching trigonometric moments is useful as a diagnostic and for starting values, but it is not
#' a reliable estimator for this family; \code{\link{fit_wsgn}} fits the full four-parameter
#' distribution by maximum likelihood.
#' @param theta numeric vector of angles in radians.
#' @param starts list of numeric starting values \code{c(mu, beta, lam)}.
#' @return a list with the estimate \code{par = c(mu, beta, lam)}, its standard errors \code{se},
#'   and the objective value \code{value}.
#' @examples
#' set.seed(11)
#' truth <- c(mu = 1, beta = 1.5, lam = 0.3)
#' th <- rwsgn(1000, truth["mu"], truth["beta"], truth["lam"])
#' fit <- fit_wsgn_trig(th)
#' round(rbind(truth = truth, estimate = fit$par, se = fit$se), 3)
#' @export
fit_wsgn_trig <- function(theta, starts = NULL) {
  m <- c(mean(cos(theta)), mean(sin(theta)), mean(cos(2 * theta)), mean(sin(2 * theta)))
  mu0 <- atan2(m[2], m[1]) %% (2 * pi)
  if (is.null(starts)) starts <- list(c(mu0, 2, 0.5), c(mu0, 1, 0.3), c(mu0, 4, 0.7))
  mom <- function(p3) {
    r <- trig_moments(2, p3[1], p3[2], p3[3])
    c(r$A[1], r$B[1], r$A[2], r$B[2])
  }
  obj <- function(p3) {
    if (p3[2] < 0.3 || p3[2] > 12 || p3[3] <= 0.02 || p3[3] >= 0.98) return(1e7)
    v <- sum((mom(p3) - m)^2)
    if (is.finite(v)) v else 1e7
  }
  fits <- lapply(starts, optim, fn = obj, method = "Nelder-Mead")
  best <- fits[[which.min(vapply(fits, function(r) r$value, 0))]]
  est <- best$par
  se <- md_sandwich_se(theta, est, mom)
  list(par = setNames(est, c("mu", "beta", "lam")),
       se = setNames(se, c("mu", "beta", "lam")), value = best$value)
}

## Minimum-distance sandwich standard errors, (G'G)^{-1} G' V G (G'G)^{-1} / n, with G the Jacobian
## of the population moments by central differences and V the covariance of the sample moments.
md_sandwich_se <- function(theta, est, mom) {
  vmat <- cov(cbind(cos(theta), sin(theta), cos(2 * theta), sin(2 * theta)))
  gmat <- vapply(1:3, function(j) {
    h <- 1e-4
    e <- rep(0, 3)
    e[j] <- h
    (mom(est + e) - mom(est - e)) / (2 * h)
  }, numeric(4))
  bmat <- solve(crossprod(gmat), t(gmat))
  sqrt(diag(bmat %*% vmat %*% t(bmat)) / length(theta))
}

#' Fit the reciprocal-rate subfamily by maximum likelihood
#'
#' Nelder-Mead over the three parameters \code{(mu, beta, lam)} of the reciprocal-rate subfamily
#' \code{sigma = sqrt(lam (1 - lam))}, a restriction of the model on which the likelihood is
#' bounded, with optional parametric-bootstrap standard errors. Kept from version 0.1.0 for backward
#' compatibility; \code{\link{fit_wsgn}} fits the full four-parameter distribution. Since version
#' 0.1.1 the wrapped density includes the Euler-Maclaurin tail, so log-likelihoods at shapes below
#' about 0.6 differ from those of version 0.1.0, which were truncated.
#' @param theta numeric vector of angles in radians.
#' @param start numeric starting value \code{c(mu, beta, lam)}; default is the trigonometric-moment
#'   estimate.
#' @param K largest number of images summed directly, as in \code{\link{dwsgn}}.
#' @param B number of parametric-bootstrap replicates for standard errors; \code{0} skips them.
#' @return a list with \code{par}, \code{logLik}, and (if \code{B > 0}) \code{se}.
#' @examples
#' set.seed(12)
#' truth <- c(mu = 1, beta = 1.5, lam = 0.3)
#' th <- rwsgn(500, truth["mu"], truth["beta"], truth["lam"])
#' fit <- fit_wsgn_mle(th)
#' round(rbind(truth = truth, estimate = fit$par), 3)
#' fit$logLik
#' \donttest{
#' ## bootstrap standard errors refit the model once per replicate
#' fit_b <- fit_wsgn_mle(th, B = 30)
#' round(rbind(estimate = fit_b$par, se = fit_b$se), 3)
#' }
#' @export
fit_wsgn_mle <- function(theta, start = NULL, K = 50, B = 0) {
  if (is.null(start)) start <- fit_wsgn_trig(theta)$par
  nll <- function(p3, th = theta) {
    if (p3[2] < 0.3 || p3[2] > 12 || p3[3] <= 0.02 || p3[3] >= 0.98) return(1e7)
    v <- .C(C_wsgn_nll, theta = as.double(th), nt = length(th), mu = as.double(p3[1]),
            sigma = as.double(sqrt(p3[3] * (1 - p3[3]))), beta = as.double(p3[2]),
            lam = as.double(p3[3]), K = as.integer(K), out = double(1))$out
    if (is.finite(v)) v else 1e7
  }
  r <- optim(start, nll, method = "Nelder-Mead")
  out <- list(par = setNames(r$par, c("mu", "beta", "lam")), logLik = -r$value)
  if (B > 0) {
    p <- out$par
    ests <- vapply(seq_len(B), function(i) {
      x <- rwsgn(length(theta), mu = p["mu"], beta = p["beta"], lam = p["lam"])
      optim(p, nll, th = x, method = "Nelder-Mead")$par
    }, numeric(3))
    out$se <- setNames(apply(ests, 1, sd), c("mu", "beta", "lam"))
  }
  out
}

## ------------------------------------------------------------ four-parameter maximum likelihood

.one_sided_gap <- log(1e12)

.grid_check <- function(resolution) {
  if (is.null(resolution)) return(invisible(NULL))
  if (length(resolution) != 1L || !is.finite(resolution) || resolution <= 0 || resolution > pi)
    stop("'resolution' must be one positive number of radians, at most pi")
  J <- round(2 * pi / resolution)
  if (J > 1e6)
    stop("'resolution' gives ", format(J, scientific = FALSE), " cells; at most 1e6 are allowed ",
         "(resolution at least 2 * pi / 1e6)")
  if (abs(J * resolution - 2 * pi) > 1e-8)
    stop("'resolution' must divide the circle into a whole number of cells, ",
         "e.g. pi / 18 for 10 degrees")
  invisible(J)
}

.fit_args <- function(theta, beta_floor, resolution, origin, symmetric, window, K) {
  if (!is.numeric(theta) || length(theta) < 3L || !all(is.finite(theta)))
    stop("'theta' must be at least three finite angles (no NA, NaN or Inf)")
  grouped <- !is.null(resolution)
  J <- .grid_check(resolution)
  if (length(origin) != 1L || !is.finite(origin)) stop("'origin' must be one finite number")
  if (grouped) {
    if (length(unique(round((theta - origin) / resolution) %% J)) < 2L)
      stop("the data occupy a single cell: the grouped likelihood has no interior maximum")
  } else if (length(unique(theta %% (2 * pi))) < 3L) {
    stop("exact data need at least three distinct angles: with fewer, the likelihood is ",
         "unbounded at every shape")
  }
  if (is.null(beta_floor)) beta_floor <- if (grouped) 1e-3 else 0.1
  if (length(beta_floor) != 1L || !is.finite(beta_floor) || beta_floor <= 0 || beta_floor >= 20)
    stop("'beta_floor' must be one number in (0, 20)")
  .check_K(K)
  list(theta = as.double(theta), n = length(theta), grouped = as.integer(grouped),
       res = as.double(if (grouped) resolution else 0), origin = as.double(origin),
       beta0 = as.double(beta_floor), sym = as.integer(isTRUE(symmetric)),
       window = as.double(window), K = as.integer(K))
}

.fit_out_names <- c("mu", "sigma", "lam", "beta", "logLik", "nfev", "at_floor", "starts",
                    "at_upper", "at_rate_bound", "log_a", "log_b")

.max_cores <- 256L
.cap_cores <- function(ncores) {
  if (ncores > .max_cores) {
    message("'ncores' reduced from ", format(ncores), " to ", .max_cores, ", the most threads used")
    ncores <- .max_cores
  }
  ncores
}

#' Fit the wrapped two-piece exponential power distribution by maximum likelihood (C back-end)
#'
#' Maximum likelihood fitting of the full four-parameter wrapped distribution
#' \code{(mu, sigma, lam, beta)}.
#' Every step runs in C: the optimiser, the starting values, the location search and the bootstrap
#' refits.
#'
#' \emph{Exact data.} With the shape free the likelihood is unbounded as \code{beta} tends to zero,
#' so it is maximised under a lower bound \code{beta >= beta_floor} (default \code{0.1}), under
#' which it is bounded and the estimator is consistent whenever the floor does not exceed the true
#' shape. The estimate is the maximiser of the likelihood over \code{beta_floor <= beta <= 20},
#' \code{|log a|, |log b| <= R}, where \code{R = max(30, 3 (1 - log beta_floor) / beta_floor)} (99
#' at the default floor). Two searches are made and the best result is returned. The first, for the
#' regular branch, runs Nelder-Mead from 27 starting values about the circular mean and 6 one-sided
#' starts (one log-rate near the bound, the location at either end of the widest gap between
#' observations); when its shape is at most one, every maximising location is an observation, so the
#' location is then searched over the stored observations within \code{window} of it, followed by a
#' free polish. The second is for the spike branch, whose likelihood grows as the shape falls and
#' whose profile in the shape can be ragged just above the floor: the same procedure, with a spike
#' start on the most-tied observation, is run at the fixed shapes \code{beta_floor} times 1, 1.25
#' and 1.6 with the location searched over every observation, each followed by a free polish (the
#' 1.25 and 1.6 shapes are skipped when they reach half the upper bound, that is when the floor is
#' 6.25 or more). In every location search each candidate observation is screened after a short optimisation of the
#' other parameters (from the current fit and from a spike as wide as the gap to its nearest
#' neighbour); all candidates are screened when there are at most 200 of them, otherwise the 15 best
#' at the current fit, the 15 with the smallest gap to a neighbour, the 15 with the smallest gap to
#' the further neighbour and every tied value. The 12 best screened candidates and every tied value
#' are refined further, and the best two fully. \code{\link{profile_wsgn}} uses the same search at
#' its fixed shapes, so the profile at the floor equals the fit's floor search. With no ties and a
#' moderate sample the spike often wins, and the fit is then on the floor or just above it:
#' \code{at_floor} reports a fit on the floor. \code{at_upper} reports a fit on the upper bound 20
#' and \code{at_rate_bound} one on the rate bound \code{R} (a degenerate fit, for instance a spike
#' or a one-sided limit, where \code{lam} may round to 0 or 1 and \code{sigma} to almost 0; the
#' log-rates in \code{log_rates} are exact). \code{at_one_sided} reports a one-sided fit that
#' stops short of the rate bound (one rate effectively infinite, \code{lam} within about
#' \code{1e-12} of 0 or 1, so the two log-rates differ by more than \code{log(1e12)}), for exact and
#' grouped data alike, judged from the log-rates; it is
#' the true maximum on the boundary, kept as such. The search is not exhaustive. On the cases checked (the
#' review cases, the periwinkles and the wind directions) it matched or beat independent searches,
#' but more distant branches cannot be excluded; compare with \code{\link{profile_wsgn}}. Without
#' standard errors a fit takes about 2.5 seconds for the 31 periwinkles, 5 to 10 seconds for 60
#' observations, 20 for 150 and 40 for 300 on a current laptop.
#'
#' Exact data need at least three distinct angles, and all angles must be finite.
#'
#' \emph{Rounded data.} When the data are recorded to a resolution (whole degrees, 5 or 10 degrees,
#' compass points), give \code{resolution} (in radians): the observations are assigned to cells of
#' that width centred at \code{origin + j * resolution}, and the grouped likelihood
#' \code{sum_j N_j log P(Theta in C_j)}, the likelihood of what was observed, is maximised. It is
#' bounded with no constraint on the shape (the default floor is \code{0.001}; the upper bound is
#' 50 and the rate bound 60), with the cell probabilities computed exactly from the incomplete-gamma
#' survival functions. Starting values are 30 fixed-seed random points, every fifth a spike on the
#' modal cell. The data must occupy at least two cells.
#'
#' \emph{Standard errors.} With \code{se = "auto"}, a continuous fit with \code{beta > 1} that is on
#' none of the bounds gets standard errors from the expected Fisher information \code{I(theta)^-1 /
#' n}, which is finite for \code{beta > 1/2}, the normal approximation being proved for a true shape
#' above one; in a simulation at \code{n = 200} the 95 per cent Wald intervals covered 87 to 92 per
#' cent, so they are somewhat optimistic at moderate \code{n}. In every other case \code{se =
#' "auto"} returns no standard errors, with a message saying why: on a bound they have no meaning,
#' and at a fitted shape at most one, or for grouped data, the bootstrap is available on request.
#' A one-sided fit (\code{at_one_sided}: \code{min(lam, 1 - lam)} below about \code{1e-12}, so
#' that one rate is effectively infinite, a maximum on the boundary that the optimiser can reach
#' short of the rate bound) counts as on a bound. Standard errors are never returned as \code{NaN}, except the \code{lam} entry of a symmetric fit, where \code{lam} is fixed:
#' if the information matrix is not finite and positive definite, \code{se = "auto"} gives none,
#' with a message. \code{se = "information"} is refused for grouped data, a fit on a bound (the floor, the upper
#' bound or the rate bound), a one-sided fit, and when the information matrix cannot be inverted. \code{se =
#' "bootstrap"} draws \code{B} samples from the fitted rates (rounded to the same cells when
#' grouped) and refits each by the same procedure; a refit takes as long as a fit, so \code{B = 200}
#' takes minutes to hours (\code{ncores} helps where the package is built with OpenMP). At a fitted
#' shape below one these are descriptive.
#'
#' \emph{Speed.} A fit is a search from many starting values, a screen of candidate locations
#' and, at the smallest shapes, a spike search, and these trials are independent. With
#' \code{ncores > 1} and a package built with OpenMP they run in parallel, and the best is chosen
#' afterwards in a fixed order (ties to the lowest index), so the fit is identical, to the last bit, for every
#' number of threads. A bootstrap spreads its refits over the threads instead, each refit running
#' on one thread. \code{ncores = 1} (the default, as CRAN requires) runs everything serially. A
#' likelihood evaluation with 128 or more observations (or cells) is itself split over the threads
#' and summed in a fixed order. On a laptop with four performance and six efficiency cores, busy
#' with other work, a fit of 1500 observations took 160 seconds on one thread, 86 on two, 57 on
#' four and 55 on eight, and one of 400 observations 36, 18, 11 and 13 seconds: about a factor of
#' 2.8 on four threads, and little more beyond. The gain stops there because some of the work
#' still runs in sequence (the first pass of the location search, the Nelder-Mead runs and the
#' final polish), and because the extra threads of such a machine run on its slower cores. A fit
#' takes minutes from about a thousand observations on. Where the package is built without OpenMP
#' the argument is accepted and has no effect.
#' @param theta numeric vector of angles in radians.
#' @param beta_floor lower bound on the shape; \code{NULL} gives \code{0.1} for exact data and
#'   \code{0.001} for grouped data.
#' @param resolution width of the recording cells in radians (for instance \code{pi / 18} for 10
#'   degrees), dividing the circle into a whole number of cells, at most \code{1e6} of them
#'   (a finer resolution is refused); \code{NULL} for exact data.
#' @param origin centre of one recording cell, in radians.
#' @param symmetric logical; fit the symmetric member \code{lam = 0.5}.
#' @param se \code{"auto"}, \code{"information"}, \code{"bootstrap"} or \code{"none"}.
#' @param B number of bootstrap replicates.
#' @param window half-width (radians) of the location search about the fitted location.
#' @param K largest number of images summed directly, as in \code{\link{dwsgn}}.
#' @param ngrid number of points of the circle in the quadrature for the information, a whole
#'   number from 100 to \code{1e8}.
#' @param ncores number of threads, a whole number of at least one (at most 256; more are reduced
#'   to 256 with a message), used when the package is built
#'   with OpenMP: for the independent trials inside a fit and, with \code{se = "bootstrap"}, for
#'   the refits. The result does not depend on it. Keep it at \code{1} in examples and tests.
#' @return an object of class \code{"wsgn_fit"}: a list with \code{par} (\code{mu}, \code{sigma},
#'   \code{lam}, \code{beta}), \code{rates} (\code{a}, \code{b}), \code{log_rates}
#'   (\code{log_a}, \code{log_b}, exact also on a degenerate fit), \code{logLik}, \code{se},
#'   \code{vcov} (information-based standard errors only), \code{boot} (the bootstrap fits, if any),
#'   \code{frac_beta_lt1} (fraction of bootstrap fits with shape below one), \code{se_method},
#'   \code{beta_floor}, \code{at_floor}, \code{at_upper}, \code{at_rate_bound},
#'   \code{at_one_sided}, \code{n},
#'   \code{grouped}, \code{resolution},
#'   \code{symmetric} and \code{nfev}.
#' @seealso \code{\link{profile_wsgn}} for the profile likelihood in the shape, and
#'   \code{\link{fit_wsgn_mle}} for the reciprocal-rate subfamily of version 0.1.0.
#' @examples
#' ## simulate from a skew member with a shape above one and recover it
#' set.seed(4)
#' th <- rwsgn(40, mu = 1, beta = 1.6, lam = 0.3, sigma = 0.4)
#' fit <- fit_wsgn(th)
#' fit
#' \donttest{
#' ## the same fit on two threads (where the package is built with OpenMP) is identical
#' identical(fit$par, fit_wsgn(th, ncores = 2)$par)
#' ## the same data recorded to whole multiples of 10 degrees: the grouped likelihood
#' h <- 10 * pi / 180
#' th10 <- (round(th / h) * h) %% (2 * pi)
#' fit10 <- fit_wsgn(th10, resolution = h, se = "none")
#' rbind(exact = coef(fit), grouped = coef(fit10))
#' ## a shape below one: the location is an observation (se = "bootstrap" would add
#' ## bootstrap standard errors, at the cost of one such fit per sample)
#' th2 <- rwsgn(60, mu = 2, beta = 0.9, lam = 0.4, sigma = 0.3)
#' fit2 <- fit_wsgn(th2, se = "none")
#' coef(fit2)
#' any(th2 %% (2 * pi) == fit2$par["mu"])
#' }
#' @export
fit_wsgn <- function(theta, beta_floor = NULL, resolution = NULL, origin = 0, symmetric = FALSE,
                     se = c("auto", "information", "bootstrap", "none"), B = 200, window = 0.4,
                     K = 50, ngrid = 20000, ncores = 1) {
  se <- match.arg(se)
  a <- .fit_args(theta, beta_floor, resolution, origin, symmetric, window, K)
  whole <- function(v, lo) {
    length(v) == 1L && is.numeric(v) && is.finite(v) && v >= lo && v == round(v)
  }
  if (!whole(B, 2) || !whole(ngrid, 100) || !whole(ncores, 1))
    stop("'B' (at least 2), 'ngrid' (at least 100) and 'ncores' must be whole numbers")
  if (ngrid > 1e8) stop("'ngrid' must be at most 1e8")
  ncores <- .cap_cores(ncores)
  if (!is.numeric(window) || length(window) != 1L || !is.finite(window) || window < 0)
    stop("'window' must be one non-negative number")
  r <- .C(C_wsgn_fit, theta = a$theta, n = a$n, grouped = a$grouped, res = a$res,
          origin = a$origin, beta0 = a$beta0, sym = a$sym, fixbeta = 0, window = a$window,
          K = a$K, nthreads = as.integer(ncores), out = double(12))$out
  r <- setNames(r, .fit_out_names)
  par <- r[c("mu", "sigma", "lam", "beta")]
  ## one-sided: min(lam, 1 - lam) below about 1e-12, i.e. one rate effectively infinite; judged
  ## from the exact log-rates, which stay valid when lam rounds to 0 or 1
  one_sided <- abs(r[["log_a"]] - r[["log_b"]]) > .one_sided_gap
  bound <- r["at_floor"] == 1 || r["at_upper"] == 1 || r["at_rate_bound"] == 1 || one_sided
  method <- se
  if (se == "auto") {
    method <- if (!a$grouped && par["beta"] > 1 && !bound) "information" else "none"
    if (method == "none")
      message("no standard errors: ",
              if (one_sided && !(r["at_floor"] == 1 || r["at_upper"] == 1 || r["at_rate_bound"] == 1))
                "the fit is one-sided (lam is numerically 0 or 1, one rate effectively infinite), on the boundary of the parameter space, where they have no meaning"
              else if (bound) "the fit is on a bound of the parameter space, where they have no meaning"
              else if (a$grouped) "the information is not used for grouped data"
              else "the fitted shape is at most one, where the normal approximation is not proved",
              if (!bound) "; se = \"bootstrap\" gives descriptive bootstrap standard errors"
              else "")
  }
  if (method == "information" && (a$grouped || bound))
    stop("information standard errors are not available for grouped data, for a fit on a bound ",
         "or for a one-sided fit (lam numerically 0 or 1); use se = \"bootstrap\"")
  out <- list(par = par,
              rates = c(a = exp(unname(r["log_a"])), b = exp(unname(r["log_b"]))),
              log_rates = c(log_a = unname(r["log_a"]), log_b = unname(r["log_b"])),
              logLik = unname(r["logLik"]), se = NULL, vcov = NULL, boot = NULL,
              frac_beta_lt1 = NULL, se_method = method, beta_floor = a$beta0,
              at_floor = as.logical(r["at_floor"]), at_upper = as.logical(r["at_upper"]),
              at_rate_bound = as.logical(r["at_rate_bound"]), at_one_sided = one_sided, n = a$n,
              grouped = as.logical(a$grouped),
              resolution = resolution, symmetric = isTRUE(symmetric), nfev = unname(r["nfev"]))
  if (method == "information") {
    s <- .C(C_wsgn_info, par = as.double(par), sym = a$sym, K = a$K, ngrid = as.integer(ngrid),
            n = a$n, info = double(16), vcov = double(16), se = double(4))
    nm <- c("mu", "sigma", "lam", "beta")
    act <- if (a$sym) c(1L, 2L, 4L) else 1:4
    I <- matrix(s$info, 4, 4)[act, act, drop = FALSE]
    pd <- all(is.finite(I)) && all(is.finite(s$se[act])) && all(s$se[act] > 0) &&
      !inherits(try(chol(I), silent = TRUE), "try-error")
    if (pd) {
      out$se <- setNames(s$se, nm)
      out$vcov <- matrix(s$vcov, 4, 4, dimnames = list(nm, nm))
    } else if (se == "information") {
      stop("the information matrix is not finite and positive definite at this fit, so there ",
           "are no information standard errors; use se = \"bootstrap\"")
    } else {
      message("no standard errors: the information matrix is not finite and positive definite ",
              "at this fit; se = \"bootstrap\" gives descriptive bootstrap standard errors")
      out$se_method <- "none"
    }
  } else if (method == "bootstrap") {
    ## drawn from the rates, which stay valid when lam rounds to 0 or 1 on a degenerate fit
    draws <- rsgn(a$n * B, mu = par["mu"], sigma = 1, beta = par["beta"],
                  lam1 = out$rates[["a"]], lam2 = out$rates[["b"]], wrap = TRUE)
    bt <- .fit_batch(draws, a, B, c(par, out$log_rates), ncores)
    out$se <- bt$se
    out$boot <- bt$est
    out$frac_beta_lt1 <- bt$frac_beta_lt1
  }
  structure(out, class = "wsgn_fit")
}

## refit each column of an n x B matrix of samples by the C procedure (also used by the parity
## harness, which passes the same samples to both front ends)
.fit_batch <- function(draws, a, B, par, ncores = 1) {
  ## par0 for the C code: (mu, sigma, lam, beta, log a, log b)
  if (length(par) == 4L)
    par <- c(par, log(par[3] / par[2]), log((1 - par[3]) / par[2]))
  s <- .C(C_wsgn_fit_batch, samples = as.double(draws), n = a$n, B = as.integer(B),
          grouped = a$grouped, res = a$res, origin = a$origin, beta0 = a$beta0, sym = a$sym,
          window = a$window, K = a$K, nthreads = as.integer(ncores), par0 = as.double(par),
          est = double(12 * B), se = double(5))
  est <- matrix(s$est, nrow = B, byrow = TRUE, dimnames = list(NULL, .fit_out_names))
  list(est = est[, c("mu", "sigma", "lam", "beta", "logLik")],
       se = setNames(s$se[1:4], c("mu", "sigma", "lam", "beta")), frac_beta_lt1 = s$se[5])
}

#' Profile log-likelihood of the wrapped distribution in the shape (C back-end)
#'
#' \code{M(beta)}, the supremum of the log-likelihood over \code{(mu, sigma, lam)} at each fixed
#' shape, for exact data or (with \code{resolution}) for the grouped likelihood. For exact data the
#' supremum tends to infinity as \code{beta} tends to zero, at the leading rate
#' \code{n log(n / (n - m)) / beta}, with \code{m} the largest number of tied observations, so the
#' profile shows where that degenerate branch overtakes the ordinary fit. At shapes up to one the
#' location is searched over the stored observations within \code{window} of the best free fit, as
#' described for \code{\link{fit_wsgn}} (every observation screened when there are at most 200);
#' the default \code{window = pi} covers all of them. The rate bound is that of
#' \code{\link{fit_wsgn}} with the same \code{beta_floor}, \code{max(30, 3 (1 - log b) / b)} with
#' \code{b} the smaller of \code{beta_floor} and the shape, so the profile at the floor equals the
#' floor search of the fit and the profile at any shape at or above the floor uses the fit's bound
#' (the bound is computed from \code{beta_floor}, not from each profiled shape).
#' With \code{ncores > 1} and a package built with OpenMP the shapes, and the trials inside each
#' fit, share the threads; the profile is identical for every number of threads.
#' @inheritParams fit_wsgn
#' @param ncores number of threads, a whole number of at least one (at most 256; more are
#'   reduced to 256 with a message), used when the package is built with OpenMP for the shapes and
#'   the trials inside each fit. The result does not depend on it.
#' @param beta numeric vector of shapes.
#' @param beta_floor the shape floor of the fit to compare with, which sets the rate bound
#'   (\code{NULL}: as in \code{\link{fit_wsgn}}).
#' @return a data frame with columns \code{beta}, \code{logLik}, \code{mu}, \code{sigma} and
#'   \code{lam}.
#' @examples
#' set.seed(22)
#' th <- rwsgn(50, mu = 1, beta = 1.2, lam = 0.35, sigma = 0.4)
#' pr <- profile_wsgn(th, beta = c(0.8, 1.2, 2), window = 0.4)
#' pr
#' plot(pr$beta, pr$logLik, type = "b", xlab = "shape", ylab = "profile log-likelihood")
#' @export
profile_wsgn <- function(theta, beta, resolution = NULL, origin = 0, symmetric = FALSE,
                         window = pi, K = 50, beta_floor = NULL, ncores = 1) {
  if (!(length(ncores) == 1L && is.numeric(ncores) && is.finite(ncores) && ncores >= 1 &&
        ncores == round(ncores)))
    stop("'ncores' must be a whole number, at least 1")
  ncores <- .cap_cores(ncores)
  if (!is.numeric(beta) || !all(is.finite(beta)) || any(beta <= 0))
    stop("'beta' must be positive and finite")
  a <- .fit_args(theta, beta_floor, resolution, origin, symmetric, window, K)
  r <- .C(C_wsgn_profile, theta = a$theta, n = a$n, grouped = a$grouped, res = a$res,
          origin = a$origin, beta0 = a$beta0, betas = as.double(beta), nb = length(beta), sym = a$sym,
          window = a$window, K = a$K, nthreads = as.integer(ncores),
          out = double(12 * length(beta)))$out
  m <- matrix(r, ncol = 12, byrow = TRUE, dimnames = list(NULL, .fit_out_names))
  data.frame(beta = beta, logLik = m[, "logLik"], mu = m[, "mu"], sigma = m[, "sigma"],
             lam = m[, "lam"])
}

#' Methods for fitted wrapped two-piece exponential power distributions
#'
#' Print, coefficients, covariance matrix and log-likelihood of a fit from \code{\link{fit_wsgn}}.
#' @param x,object a \code{"wsgn_fit"} object.
#' @param digits number of significant digits to print.
#' @param ... unused.
#' @return \code{print} returns \code{x} invisibly; \code{coef} the estimates
#'   \code{(mu, sigma, lam, beta)}; \code{vcov} the information-based covariance matrix
#'   (\code{NULL} for bootstrap standard errors, whose replicates are in \code{object$boot});
#'   \code{logLik} an object of class \code{"logLik"} (the grouped log-likelihood for grouped data).
#' @examples
#' set.seed(7)
#' fit <- fit_wsgn(rwsgn(30, mu = 1, beta = 2, lam = 0.4, sigma = 0.5))
#' print(fit)
#' coef(fit)
#' vcov(fit)
#' logLik(fit)
#' AIC(fit)
#' @name wsgn_fit-methods
NULL

#' @rdname wsgn_fit-methods
#' @export
print.wsgn_fit <- function(x, digits = 4, ...) {
  cat("Wrapped two-piece exponential power distribution, maximum likelihood",
      if (x$grouped) "(grouped)" else "", "\n")
  cat("n =", x$n, " log-likelihood =", format(x$logLik, digits = digits + 3),
      " shape floor =", format(x$beta_floor), if (x$at_floor) "(fit ON the floor)" else "",
      if (isTRUE(x$at_upper)) "(fit ON the upper bound 20)" else "",
      if (isTRUE(x$at_rate_bound)) "(fit ON the rate bound: degenerate)" else "",
      if (isTRUE(x$at_one_sided)) "(one-sided fit: lam numerically 0 or 1)" else "", "\n")
  tab <- rbind(estimate = x$par, se = if (is.null(x$se)) NA else x$se)
  print(signif(tab, digits))
  cat("standard errors:", x$se_method, "\n")
  invisible(x)
}

#' @rdname wsgn_fit-methods
#' @export
coef.wsgn_fit <- function(object, ...) object$par

#' @rdname wsgn_fit-methods
#' @export
vcov.wsgn_fit <- function(object, ...) object$vcov

#' @rdname wsgn_fit-methods
#' @export
logLik.wsgn_fit <- function(object, ...) {
  structure(object$logLik, df = if (object$symmetric) 3L else 4L, nobs = object$n,
            class = "logLik")
}

#' Grouped log-likelihood of rounded circular data (C back-end)
#'
#' The likelihood of what is observed when angles are recorded to a resolution: the observations
#' are assigned to cells of width \code{resolution} centred at \code{origin + j * resolution}, and
#' the grouped log-likelihood is \code{sum_j N_j log P(Theta in C_j)}, with the cell probabilities
#' computed exactly as in \code{\link{pwsgn}}. It is at most zero and bounded for every parameter
#' value, and \code{\link{fit_wsgn}} maximises it when \code{resolution} is given.
#' @inheritParams fit_wsgn
#' @param mu,sigma,lam,beta parameters of the wrapped distribution, as in \code{\link{dwsgn}}.
#' @return the grouped log-likelihood.
#' @examples
#' ## ants' directions recorded to the nearest 10 degrees (Fisher, 1993, data set B.7)
#' if (requireNamespace("circular", quietly = TRUE)) {
#'   e <- new.env()
#'   utils::data("fisherB7", package = "circular", envir = e)
#'   ants <- e$fisherB7 * pi / 180
#'   grouped_loglik(ants, pi / 18, mu = pi, sigma = 0.1, lam = 0.5, beta = 0.6)
#' }
#' ## the grouped likelihood of a single cell holding all the data is zero at a spike
#' grouped_loglik(rep(1, 20), pi / 18, mu = 1, sigma = 1e-3, lam = 0.5, beta = 2, origin = 1)
#' @export
grouped_loglik <- function(theta, resolution, mu, sigma, lam, beta, origin = 0, K = 50) {
  bad_res <- missing(resolution) || is.null(resolution) || !is.numeric(resolution) ||
    length(resolution) != 1L
  if (bad_res)
    stop("'resolution' must be given, as one number of radians")
  if (!is.numeric(theta) || !all(is.finite(theta)))
    stop("'theta' must be finite angles (no NA, NaN or Inf)")
  if (length(origin) != 1L || !is.finite(origin)) stop("'origin' must be one finite number")
  .grid_check(resolution)
  .check_wrapped(sigma, beta, lam, K)
  .C(C_wsgn_grouped_ll, theta = as.double(theta), n = length(theta), res = as.double(resolution),
     origin = as.double(origin), mu = as.double(mu), sigma = as.double(sigma),
     beta = as.double(beta), lam = as.double(lam), K = as.integer(K), out = double(1))$out
}
