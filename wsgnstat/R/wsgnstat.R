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
#' x <- seq(-4, 4, length.out = 9)
#' round(dsgn(x, 0, 1, 2, 0.5, 0.5), 5)
#' @export
dsgn <- function(x, mu = 0, sigma = 1, beta = 2, lam1 = 0.5, lam2 = 0.5) {
  .C(C_wsgn_density, x = as.double(x), nx = length(x), mu = as.double(mu),
     sigma = as.double(sigma), beta = as.double(beta), lam1 = as.double(lam1),
     lam2 = as.double(lam2), out = double(length(x)))$out
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
#' @param K truncation of the wrapping sum (\code{abs(k) <= K}).
#' @return the wrapped probability density function evaluated at \code{theta}.
#' @examples
#' ## at beta = 2 with equal weights the circular law IS the wrapped normal
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

#' Characteristic function of the skew generalised normal distribution (C back-end)
#'
#' Evaluated through the absolutely convergent cosine and sine integrals underlying the Fox
#' H-function representation, hence valid for every \code{beta > 0}.
#' @param t numeric vector of arguments.
#' @inheritParams dsgn
#' @return a complex vector of the same length as \code{t}.
#' @examples
#' ## reached through absolutely convergent integrals, so every beta > 0 is available
#' cfsgn(0, 0, 1, 2, 0.5, 0.5)
#' Mod(cfsgn(c(0.5, 1.5, 3), 0, 1, 1.3, 0.8, 0.3))
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
#' ## the same reduction in the trigonometric moments: exp(-p^2 / 4)
#' tm <- trig_moments(4, mu = 0, beta = 2, lam = 0.5)
#' rbind(engine = tm$A, closed_form = exp(-(1:4)^2 / 4))
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
#' @export
rsgn <- function(n, mu = 0, sigma = 1, beta = 2, lam1 = 0.5, lam2 = 0.5, wrap = FALSE) {
  u <- rgamma(n, shape = 1 + 1/beta)^(1/beta)
  x <- runif(n, min = mu - (sigma/lam2)*u, max = mu + (sigma/lam1)*u)
  if (wrap) x %% (2*pi) else x
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
#' @export
fit_wsgn_trig <- function(theta, starts = NULL) {
  m <- c(mean(cos(theta)), mean(sin(theta)), mean(cos(2*theta)), mean(sin(2*theta)))
  mu0 <- atan2(m[2], m[1]) %% (2*pi)
  if (is.null(starts)) starts <- list(c(mu0, 2, 0.5), c(mu0, 1, 0.3), c(mu0, 4, 0.7))
  mom <- function(p3) { r <- trig_moments(2, p3[1], p3[2], p3[3]); c(r$A[1], r$B[1], r$A[2], r$B[2]) }
  obj <- function(p3) {
    if (p3[2] < 0.3 || p3[2] > 12 || p3[3] <= 0.02 || p3[3] >= 0.98) return(1e7)
    v <- sum((mom(p3) - m)^2); if (is.finite(v)) v else 1e7
  }
  best <- NULL
  for (s in starts) { r <- optim(s, obj, method = "Nelder-Mead")
    if (is.null(best) || r$value < best$value) best <- r }
  est <- best$par
  V <- cov(cbind(cos(theta), sin(theta), cos(2*theta), sin(2*theta)))
  G <- sapply(1:3, function(j) { h <- 1e-4; e <- rep(0, 3); e[j] <- h
    (mom(est + e) - mom(est - e))/(2*h) })
  Bm <- solve(crossprod(G), t(G))
  se <- sqrt(diag(Bm %*% V %*% t(Bm))/length(theta))
  list(par = setNames(est, c("mu", "beta", "lam")),
       se = setNames(se, c("mu", "beta", "lam")), value = best$value)
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
#' @export
fit_wsgn_mle <- function(theta, start = NULL, K = 50, B = 0) {
  if (is.null(start)) start <- fit_wsgn_trig(theta)$par
  nll <- function(p3) {
    if (p3[2] < 0.3 || p3[2] > 12 || p3[3] <= 0.02 || p3[3] >= 0.98) return(1e7)
    v <- .C(C_wsgn_nll, theta = as.double(theta), nt = length(theta), mu = as.double(p3[1]),
            beta = as.double(p3[2]), lam = as.double(p3[3]), K = as.integer(K),
            out = double(1))$out
    if (is.finite(v)) v else 1e7
  }
  r <- optim(start, nll, method = "Nelder-Mead")
  out <- list(par = setNames(r$par, c("mu", "beta", "lam")), logLik = -r$value)
  if (B > 0) {
    p <- out$par; sig <- sqrt(p["lam"]*(1 - p["lam"]))
    ests <- replicate(B, {
      x <- rsgn(length(theta), mu = p["mu"], sigma = sig, beta = p["beta"],
                lam1 = p["lam"], lam2 = 1 - p["lam"], wrap = TRUE)
      optim(p, nll, method = "Nelder-Mead")$par })
    out$se <- setNames(apply(ests, 1, sd), c("mu", "beta", "lam"))
  }
  out
}
