/* wsgn.h : shared pure-C back-end for the two-piece exponential power distribution on the line and
 * its wrapped version on the circle, treated as a uniform scale mixture with generalised gamma
 * mixing (the package keeps the authors' earlier name, wrapped skew generalised normal). Pure libm;
 * bound identically from R (.C) and Python (ctypes).
 *
 * Wrapped family: location mu, scale sigma > 0, skewness 0 < lam < 1, shape beta > 0; right rate
 * a = lam / sigma, left rate b = (1 - lam) / sigma, P(X < mu) = lam. The line functions take the
 * general two-piece form with right rate lam1 / sigma and left rate lam2 / sigma. */
#ifndef WSGN_H
#define WSGN_H
#ifdef __cplusplus
extern "C" {
#endif

/* Density of the line distribution at x (length nx), into out. */
void wsgn_density(const double *x, const int *nx, const double *mu, const double *sigma,
                  const double *beta, const double *lam1, const double *lam2, double *out);

/* Distribution function of the line distribution, through the regularised incomplete gamma
 * function. */
void wsgn_pdist(const double *x, const int *nx, const double *mu, const double *sigma,
                const double *beta, const double *lam1, const double *lam2, double *out);

/* Characteristic function of the line distribution at t (length nt) into re, im. method 0: the
 * integral G(w) on a rotated path, trapezoid rule in log-radius, tolerance tol relative to the L1
 * norm of each half's integrand (so the error of phi is absolute, of order tol min(1, 1/|w|));
 * method 1: the ascending (beta > 1) or descending (beta < 1) series, NaN unless its estimated
 * rounding error is below tol relative to each half. */
void wsgn_cf(const double *t, const int *nt, const double *mu, const double *sigma,
             const double *beta, const double *lam1, const double *lam2, const double *tol,
             const int *method, double *re, double *im);

/* Trigonometric moments A_p, B_p (real and imaginary parts of phi_p), p = 1..P, of the wrapped
 * distribution. */
void wsgn_trigmom(const int *P, const double *mu, const double *sigma, const double *beta,
                  const double *lam, const double *tol, double *A, double *B);

/* Decay constant C and exponent 1 + beta of |phi_p| ~ C p^-(1 + beta) as p grows, into
 * out[0], out[1]. */
void wsgn_decay(const double *sigma, const double *beta, const double *lam, double *out);

/* Wrapped density (or its logarithm when give_log) at theta: images summed directly on each side
 * of mu until the Euler-Maclaurin tail (exact incomplete-gamma integral, Bernoulli terms through
 * B_12, in log space) is accurate, trying after 8, 16, 32 and K images and after each further one. */
void wsgn_wdensity(const double *theta, const int *nt, const double *mu, const double *sigma,
                   const double *beta, const double *lam, const int *K, const int *give_log,
                   double *out);

/* Negative log-likelihood of a wrapped sample, into *out. */
void wsgn_nll(const double *theta, const int *nt, const double *mu, const double *sigma,
              const double *beta, const double *lam, const int *K, double *out);

/* Wrapped distribution function on [0, 2 pi): probability of the arc [0, theta]. */
void wsgn_wcdf(const double *theta, const int *nt, const double *mu, const double *sigma,
               const double *beta, const double *lam, const int *K, double *out);

/* Wrapped quantile function, by bisection on wsgn_wcdf, geometric while the bracket is wide. */
void wsgn_wquant(const double *p, const int *np, const double *mu, const double *sigma,
                 const double *beta, const double *lam, const int *K, double *out);

/* Maximum likelihood fit of the four-parameter wrapped distribution (grouped = 0: exact
 * observations under the shape floor beta0; grouped = 1: counts on the regular grid of cells of
 * width res centred at origin + j res); sym = 1 fixes lam = 1/2; fixbeta > 0 fixes the shape.
 * For exact data the search covers [beta0, 20] including the floor branch; grouped [beta0, 50].
 * out[0..9] = mu, sigma, lam, beta, log-likelihood, evaluations, at-floor flag, starts,
 * at-upper-bound flag, at-rate-bound flag. When the library is built with OpenMP, the independent
 * trials inside the fit run on *nthreads threads (1 = serial) with a result that does not depend on
 * the number of threads. */
void wsgn_fit(const double *theta, const int *n, const int *grouped, const double *res,
              const double *origin, const double *beta0, const int *sym, const double *fixbeta,
              const double *window, const int *K, const int *nthreads, double *out);

/* Profile log-likelihood M(beta) at each of betas (length nb), out 12 x nb as for wsgn_fit;
 * beta0 is the shape floor of the fit to compare with, which sets the rate bound. Grouped data
 * with more than 1e6 cells give NaN. *nthreads as for wsgn_fit; the shapes are spread over the
 * threads as well. */
void wsgn_profile(const double *theta, const int *n, const int *grouped, const double *res,
                  const double *origin, const double *beta0, const double *betas, const int *nb,
                  const int *sym, const double *window, const int *K, const int *nthreads,
                  double *out);

/* Grouped log-likelihood of theta binned on the grid, into *out. */
void wsgn_grouped_ll(const double *theta, const int *n, const double *res, const double *origin,
                     const double *mu, const double *sigma, const double *beta, const double *lam,
                     const int *K, double *out);

/* Refit each column of samples (n x B) by the same procedure; est 10 x B; se[0..3] bootstrap
 * standard errors of (mu, sigma, lam, beta), se[4] fraction of refits with beta < 1. */
void wsgn_fit_batch(const double *samples, const int *n, const int *B, const int *grouped,
                    const double *res, const double *origin, const double *beta0, const int *sym,
                    const double *window, const int *K, const int *nthreads, const double *par0,
                    double *est, double *se);

/* Expected Fisher information per observation in (mu, sigma, lam, beta) (4 x 4, column-major),
 * vcov = I^-1 / n and standard errors se. */
void wsgn_info(const double *par, const int *sym, const int *K, const int *ngrid, const int *n,
               double *info, double *vcov, double *se);

#ifdef __cplusplus
}
#endif
#endif
