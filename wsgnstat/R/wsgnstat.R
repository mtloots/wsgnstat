## wsgnstat: the skew generalised normal distribution on the line and its wrapped version on the
## circle (Loots-Bekker), on one shared pure-C back-end (compiled at install), also bound from
## Python (wsgnstat). The characteristic function is evaluated through absolutely convergent
## integrals (the Fox H route), so densities, trigonometric moments and both estimators are
## available for every shape parameter beta > 0, not only where the Fox-Wright series converges.
#' @useDynLib wsgnstat, .registration = TRUE, .fixes = "C_"
#' @importFrom stats optim rgamma runif sd cov setNames
#' @keywords internal
"_PACKAGE"

#' Density of the skew generalised normal distribution on the line (C back-end)
#'
#' The scale mixture of a uniform and a generalised gamma variable, with two-piece exponential-power
#' form: left rate \code{lam2}, right rate \code{lam1}.
#' @param x numeric vector of quantiles.
#' @param mu location parameter.
#' @param sigma scale parameter (> 0).
#' @param beta shape (tail) parameter (> 0).
#' @param lam1 right-hand rate parameter (> 0).
#' @param lam2 left-hand rate parameter (> 0).
#' @return the probability density function evaluated at \code{x}.
#' @examples
#' ## a right-skew, heavy-tailed member and the normal member it generalises
#' x <- seq(-4, 6, length.out = 401)
#' f_skew <- dsgn(x, mu = 0, sigma = 1, beta = 1.3, lam1 = 0.3, lam2 = 0.8)
#' f_norm <- dsgn(x, mu = 0, sigma = 0.5, beta = 2, lam1 = 0.5, lam2 = 0.5)
#' plot(x, f_skew, type = "l", ylab = "density")
#' lines(x, f_norm, lty = 2)
#' ## at beta = 2 with lam1 = lam2 = 0.5 and sigma = 0.5 the law is N(0, 1/2)
#' max(abs(f_norm - dnorm(x, 0, sqrt(1 / 2))))
#' ## the density integrates to one (Riemann sum on a wide grid)
#' xx <- seq(-40, 40, length.out = 80001)
#' sum(dsgn(xx, 0, 1, 1.3, 0.3, 0.8)) * (xx[2] - xx[1])
#' @export
dsgn <- function(x, mu = 0, sigma = 1, beta = 2, lam1 = 0.5, lam2 = 0.5) {
  .C(C_wsgn_density, x = as.double(x), nx = length(x), mu = as.double(mu),
     sigma = as.double(sigma), beta = as.double(beta), lam1 = as.double(lam1),
     lam2 = as.double(lam2), out = double(length(x)))$out
}

#' Distribution function of the skew generalised normal distribution on the line (C back-end)
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
  .C(C_wsgn_pdist, x = as.double(q), nx = length(q), mu = as.double(mu),
     sigma = as.double(sigma), beta = as.double(beta), lam1 = as.double(lam1),
     lam2 = as.double(lam2), out = double(length(q)), NAOK = TRUE)$out
}

#' Density of the wrapped skew generalised normal distribution (C back-end)
#'
#' The identified parameterisation is used: \code{lam1 = lam}, \code{lam2 = 1 - lam} and
#' \code{sigma = sqrt(lam * (1 - lam))}, so the free parameters are the location \code{mu},
#' the tail parameter \code{beta} and the skewness \code{lam}.
#' @param theta numeric vector of angles in radians on \code{[0, 2*pi)}.
#' @param mu location (junction) parameter; for a skew law this is not the mean direction.
#' @param beta shape (tail) parameter (> 0).
#' @param lam skewness parameter in (0, 1); \code{0.5} is symmetric.
#' @param K truncation of the wrapping sum (\code{abs(k) <= K}). The default \code{50} captures the
#'   whole mass to double precision for \code{beta >= 0.8}; for very heavy tails (\code{beta} near
#'   \code{0.3}) it does not, and \code{pwsgn(2 * pi, mu, beta, lam, K)} reports the mass that the
#'   truncated sum captures, so \code{K} can be raised until that is one.
#' @return the wrapped probability density function evaluated at \code{theta}.
#' @examples
#' ## a skew member drawn on the circle as a radial curve about the unit circle
#' th <- seq(0, 2 * pi, length.out = 721)
#' d <- dwsgn(th, mu = pi / 2, beta = 1.5, lam = 0.3)
#' plot(cos(th), sin(th), type = "l", asp = 1, xlim = c(-2, 2), ylim = c(-2, 2),
#'      xlab = "", ylab = "", lty = 3)
#' lines((1 + d) * cos(th), (1 + d) * sin(th))
#' ## the density integrates to one over the circle
#' sum(d[-1]) * (th[2] - th[1])
#' ## at beta = 2 with lam = 0.5 the circular law IS the wrapped normal
#' ## with sigma = 1/sqrt(2), which anchors the wrapping machinery exactly
#' tt <- c(0.5, 2.0, 4.5)
#' wn <- vapply(tt, function(t) sum(dnorm(t + 2 * pi * (-80:80), 0, 1 / sqrt(2))), 0)
#' rbind(wsgnstat = dwsgn(tt, 0, 2, 0.5), wrapped_normal = wn)
#' @export
dwsgn <- function(theta, mu = 0, beta = 2, lam = 0.5, K = 50) {
  .C(C_wsgn_wdensity, theta = as.double(theta), nt = length(theta), mu = as.double(mu),
     beta = as.double(beta), lam = as.double(lam), K = as.integer(K),
     out = double(length(theta)))$out
}

#' Distribution and quantile functions of the wrapped skew generalised normal distribution
#'
#' The distribution function is measured anticlockwise from angle \code{0}:
#' \code{pwsgn(theta)} is the probability of the arc \code{[0, theta]}, computed as the sum over
#' turns of the line law's mass on \code{[2 pi k, theta + 2 pi k]} through the regularised
#' incomplete gamma function, so no numerical quadrature is involved. The quantile function inverts
#' it by bisection. The identified parameterisation of \code{\link{dwsgn}} is used.
#' @param q numeric vector of angles in radians.
#' @param p numeric vector of probabilities.
#' @inheritParams dwsgn
#' @return \code{pwsgn} gives probabilities in \code{[0, 1]}: \code{0} for \code{q <= 0} and, for
#'   \code{q >= 2 * pi}, the mass captured by the wrapping sum truncated at \code{K} (one to double
#'   precision unless the tails are very heavy; see \code{\link{dwsgn}}). \code{qwsgn} gives angles
#'   in \code{[0, 2 * pi]}, with a warning for any \code{p} above the captured mass.
#' @examples
#' ## the distribution function is the integral of the density
#' pwsgn(2, mu = 1, beta = 1.5, lam = 0.3)
#' integrate(dwsgn, 0, 2, mu = 1, beta = 1.5, lam = 0.3)$value
#' ## the quantile function inverts it
#' p <- c(0.05, 0.25, 0.5, 0.75, 0.95)
#' q <- qwsgn(p, mu = 1, beta = 1.5, lam = 0.3)
#' rbind(p, pwsgn(q, mu = 1, beta = 1.5, lam = 0.3))
#' ## the curve rises from 0 to 1 around the circle
#' th <- seq(0, 2 * pi, length.out = 361)
#' plot(th, pwsgn(th, mu = 1, beta = 1.5, lam = 0.3), type = "l",
#'      xlab = "angle (radians)", ylab = "distribution function")
#' @name pwsgn
NULL

#' @rdname pwsgn
#' @export
pwsgn <- function(q, mu = 0, beta = 2, lam = 0.5, K = 50) {
  .C(C_wsgn_wcdf, theta = as.double(q), nt = length(q), mu = as.double(mu),
     beta = as.double(beta), lam = as.double(lam), K = as.integer(K),
     out = double(length(q)))$out
}

#' @rdname pwsgn
#' @export
qwsgn <- function(p, mu = 0, beta = 2, lam = 0.5, K = 50) {
  tot <- pwsgn(2 * pi, mu, beta, lam, K)
  if (any(p > tot & p < 1, na.rm = TRUE))
    warning("the wrapping sum truncated at K = ", K, " captures mass ", format(tot),
            "; raise K for probabilities above it")
  .C(C_wsgn_wquant, p = as.double(p), np = length(p), mu = as.double(mu),
     beta = as.double(beta), lam = as.double(lam), K = as.integer(K),
     out = double(length(p)))$out
}

#' Characteristic function of the skew generalised normal distribution (C back-end)
#'
#' Evaluated through the absolutely convergent cosine and sine integrals underlying the Fox
#' H-function representation, hence valid for every \code{beta > 0}.
#' @param t numeric vector of arguments.
#' @inheritParams dsgn
#' @return a complex vector of the same length as \code{t}.
#' @examples
#' ## reached through absolutely convergent integrals, so every beta > 0 is available,
#' ## including beta <= 1, where the Fox-Wright series does not converge for all t
#' cfsgn(0, 0, 1, 2, 0.5, 0.5)
#' ## the normal member: exp(-t^2 / 4) for N(0, 1/2)
#' t <- c(0.5, 1.5, 3)
#' rbind(engine = Re(cfsgn(t, 0, 0.5, 2, 0.5, 0.5)), closed_form = exp(-t^2 / 4))
#' ## a skew member with beta = 0.8, against the Fourier transform of the density
#' x <- seq(-60, 60, length.out = 120001)
#' f <- dsgn(x, 0, 1, 0.8, 0.8, 0.3)
#' c(engine = cfsgn(1.1, 0, 1, 0.8, 0.8, 0.3),
#'   quadrature = sum(exp(1i * 1.1 * x) * f) * (x[2] - x[1]))
#' @export
cfsgn <- function(t, mu = 0, sigma = 1, beta = 2, lam1 = 0.5, lam2 = 0.5) {
  r <- .C(C_wsgn_cf, t = as.double(t), nt = length(t), mu = as.double(mu),
          sigma = as.double(sigma), beta = as.double(beta), lam1 = as.double(lam1),
          lam2 = as.double(lam2), re = double(length(t)), im = double(length(t)))
  complex(real = r$re, imaginary = r$im)
}

#' Population trigonometric moments of the wrapped family (C back-end)
#'
#' @param P highest order p.
#' @inheritParams dwsgn
#' @return a list with components \code{A} and \code{B}, the cosine and sine moments for
#'   \code{p = 1, ..., P}.
#' @examples
#' ## the reduction to the wrapped normal in the trigonometric moments: exp(-p^2 / 4)
#' tm <- trig_moments(4, mu = 0, beta = 2, lam = 0.5)
#' rbind(engine = tm$A, closed_form = exp(-(1:4)^2 / 4))
#' ## a skew member: the first moment against a large simulated sample
#' tm1 <- trig_moments(1, mu = 1, beta = 1.5, lam = 0.3)
#' set.seed(1)
#' th <- rwsgn(1e5, mu = 1, beta = 1.5, lam = 0.3)
#' rbind(population = c(tm1$A, tm1$B), sample = c(mean(cos(th)), mean(sin(th))))
#' @export
trig_moments <- function(P = 2, mu = 0, beta = 2, lam = 0.5) {
  r <- .C(C_wsgn_trigmom, P = as.integer(P), mu = as.double(mu), beta = as.double(beta),
          lam = as.double(lam), A = double(P), B = double(P))
  list(A = r$A, B = r$B)
}

#' Random generation from the (wrapped) skew generalised normal distribution
#'
#' Uses the stochastic representation: \code{U^beta} is Gamma with shape \code{1 + 1/beta}, a
#' uniform variate is drawn conditionally, and the product is wrapped when \code{wrap = TRUE}.
#' @param n sample size.
#' @inheritParams dsgn
#' @param wrap logical; reduce modulo \code{2*pi}.
#' @return a numeric vector of draws.
#' @seealso \code{\link{rwsgn}} for draws from the wrapped law in its identified
#'   parameterisation.
#' @examples
#' set.seed(2)
#' x <- rsgn(5000, mu = 0, sigma = 1, beta = 1.3, lam1 = 0.3, lam2 = 0.8)
#' hist(x, breaks = 60, freq = FALSE, main = "")
#' curve(dsgn(x, 0, 1, 1.3, 0.3, 0.8), add = TRUE)
#' ## mass left of the location: lam1 / (lam1 + lam2) = 0.2727
#' mean(x < 0)
#' @export
rsgn <- function(n, mu = 0, sigma = 1, beta = 2, lam1 = 0.5, lam2 = 0.5, wrap = FALSE) {
  u <- rgamma(n, shape = 1 + 1 / beta)^(1 / beta)
  x <- runif(n, min = mu - (sigma / lam2) * u, max = mu + (sigma / lam1) * u)
  if (wrap) x %% (2 * pi) else x
}

#' Random generation from the wrapped skew generalised normal distribution
#'
#' Draws from the circular law in the identified parameterisation of \code{\link{dwsgn}}, that is
#' \code{rsgn(n, mu, sigma = sqrt(lam * (1 - lam)), beta, lam1 = lam, lam2 = 1 - lam, wrap = TRUE)}.
#' The estimators \code{\link{fit_wsgn_trig}} and \code{\link{fit_wsgn_mle}} estimate exactly
#' these \code{(mu, beta, lam)}.
#' @inheritParams dwsgn
#' @param n number of draws.
#' @return a numeric vector of angles in \code{[0, 2 * pi)}.
#' @examples
#' set.seed(3)
#' th <- rwsgn(2000, mu = 1, beta = 1.5, lam = 0.3)
#' hist(th, breaks = seq(0, 2 * pi, length.out = 49), freq = FALSE, main = "",
#'      xlab = "angle (radians)")
#' curve(dwsgn(x, mu = 1, beta = 1.5, lam = 0.3), 0, 2 * pi, add = TRUE)
#' @export
rwsgn <- function(n, mu = 0, beta = 2, lam = 0.5) {
  rsgn(n, mu = mu, sigma = sqrt(lam * (1 - lam)), beta = beta, lam1 = lam, lam2 = 1 - lam,
       wrap = TRUE)
}

#' Fit the wrapped skew generalised normal by trigonometric moments (C back-end)
#'
#' Minimum-distance matching of the first two sample trigonometric moments to their population
#' counterparts, evaluated in C through the integral characteristic function, hence available for
#' every \code{beta > 0}. Standard errors come from the minimum-distance sandwich.
#' @param theta numeric vector of angles in radians.
#' @param starts list of numeric starting values \code{c(mu, beta, lam)}.
#' @return a list with the estimate \code{par = c(mu, beta, lam)}, its standard errors \code{se},
#'   and the objective value \code{value}.
#' @examples
#' ## simulate from a known member and recover it
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
  best <- NULL
  for (s in starts) {
    r <- optim(s, obj, method = "Nelder-Mead")
    if (is.null(best) || r$value < best$value) best <- r
  }
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

#' Fit the wrapped skew generalised normal by maximum likelihood (C back-end)
#'
#' Nelder-Mead over the identified parameters, with the wrapped log-likelihood evaluated in C.
#' Standard errors, when requested, come from a parametric bootstrap, the observed information
#' being unstable near the interior cusp the density develops as \code{beta} approaches one.
#' @param theta numeric vector of angles in radians.
#' @param start numeric starting value \code{c(mu, beta, lam)}; default is the trigonometric-moment
#'   estimate.
#' @param K truncation of the wrapping sum.
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
#' ## bootstrap standard errors refit the model once per replicate (a few seconds for B = 30)
#' fit_b <- fit_wsgn_mle(th, B = 30)
#' round(rbind(estimate = fit_b$par, se = fit_b$se), 3)
#' }
#' @export
fit_wsgn_mle <- function(theta, start = NULL, K = 50, B = 0) {
  if (is.null(start)) start <- fit_wsgn_trig(theta)$par
  nll <- function(p3, th = theta) {
    if (p3[2] < 0.3 || p3[2] > 12 || p3[3] <= 0.02 || p3[3] >= 0.98) return(1e7)
    v <- .C(C_wsgn_nll, theta = as.double(th), nt = length(th), mu = as.double(p3[1]),
            beta = as.double(p3[2]), lam = as.double(p3[3]), K = as.integer(K),
            out = double(1))$out
    if (is.finite(v)) v else 1e7
  }
  r <- optim(start, nll, method = "Nelder-Mead")
  out <- list(par = setNames(r$par, c("mu", "beta", "lam")), logLik = -r$value)
  if (B > 0) {
    p <- out$par
    ests <- replicate(B, {
      x <- rwsgn(length(theta), mu = p["mu"], beta = p["beta"], lam = p["lam"])
      optim(p, nll, th = x, method = "Nelder-Mead")$par
    })
    out$se <- setNames(apply(ests, 1, sd), c("mu", "beta", "lam"))
  }
  out
}
