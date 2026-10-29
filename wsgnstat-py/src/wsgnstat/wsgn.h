/* wsgn.h : shared pure-C back-end for the skew generalised normal distribution on the line and
 * its wrapped version on the circle (Loots-Bekker). Pure libm; bound identically from R and Python.
 * The line law is the scale mixture of a uniform and a generalised gamma variable; the wrapped law
 * uses the identified parameterisation lambda1 = lambda, lambda2 = 1 - lambda,
 * sigma = sqrt(lambda(1-lambda)). The characteristic function is evaluated through the absolutely
 * convergent integrals g_c(w) = int cos(wu) exp(-u^beta) du and g_s (the Fox H route), valid for
 * every beta > 0. */
#ifndef WSGN_H
#define WSGN_H
#ifdef __cplusplus
extern "C" {
#endif

/* Density of the line law at x (length nx), parameters mu, sigma > 0, beta > 0,
 * lam1 > 0, lam2 > 0, into out. */
void wsgn_density(const double *x, const int *nx, const double *mu, const double *sigma,
                  const double *beta, const double *lam1, const double *lam2, double *out);

/* Wrapped density on [0, 2*pi) at theta (length nt) under the identified parameterisation
 * (mu, beta, lam), truncating the wrapping sum at |k| <= K, into out. */
void wsgn_wdensity(const double *theta, const int *nt, const double *mu, const double *beta,
                   const double *lam, const int *K, double *out);

/* Characteristic function of the line law at t (length nt) into re, im, through the integral
 * representation; valid for every beta > 0. */
void wsgn_cf(const double *t, const int *nt, const double *mu, const double *sigma,
             const double *beta, const double *lam1, const double *lam2,
             double *re, double *im);

/* Population trigonometric moments A_p, B_p, p = 1..P, of the wrapped identified family. */
void wsgn_trigmom(const int *P, const double *mu, const double *beta, const double *lam,
                  double *A, double *B);

/* Negative log-likelihood of a wrapped sample theta (length nt) under the identified
 * parameterisation, wrapping sum truncated at |k| <= K, into *out. */
void wsgn_nll(const double *theta, const int *nt, const double *mu, const double *beta,
              const double *lam, const int *K, double *out);

/* Distribution function of the line law at x (length nx), through the regularised incomplete
 * gamma function, into out. */
void wsgn_pdist(const double *x, const int *nx, const double *mu, const double *sigma,
                const double *beta, const double *lam1, const double *lam2, double *out);

/* Wrapped distribution function on [0, 2*pi), measured anticlockwise from angle 0, under the
 * identified parameterisation, wrapping sum truncated at |k| <= K, into out. */
void wsgn_wcdf(const double *theta, const int *nt, const double *mu, const double *beta,
               const double *lam, const int *K, double *out);

/* Wrapped quantile function: the angle in [0, 2*pi) at which wsgn_wcdf reaches p, by bisection. */
void wsgn_wquant(const double *p, const int *np, const double *mu, const double *beta,
                 const double *lam, const int *K, double *out);

#ifdef __cplusplus
}
#endif
#endif
