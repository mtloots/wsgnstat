/* wsgn.c : shared pure-C back-end for the (wrapped) skew generalised normal distribution.
 * See wsgn.h for the interface. Pure libm so the same source compiles under R (via .C) and
 * Python (via ctypes) with byte-identical results. */
#include <math.h>
#include "wsgn.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* normalising constant of the line density */
static double norm_const(double sigma, double beta, double lam1, double lam2)
{
    return 1.0 / (sigma * (1.0 / lam1 + 1.0 / lam2) * tgamma((beta + 1.0) / beta));
}

static double line_density(double x, double mu, double sigma, double beta,
                           double lam1, double lam2, double c)
{
    double d = x - mu;
    double s = (d < 0.0) ? lam2 : lam1;
    return c * exp(-pow(s * fabs(d) / sigma, beta));
}

void wsgn_density(const double *x, const int *nx, const double *mu, const double *sigma,
                  const double *beta, const double *lam1, const double *lam2, double *out)
{
    double c = norm_const(*sigma, *beta, *lam1, *lam2);
    for (int i = 0; i < *nx; i++)
        out[i] = line_density(x[i], *mu, *sigma, *beta, *lam1, *lam2, c);
}

void wsgn_wdensity(const double *theta, const int *nt, const double *mu, const double *beta,
                   const double *lam, const int *K, double *out)
{
    double sigma = sqrt(*lam * (1.0 - *lam));
    double c = norm_const(sigma, *beta, *lam, 1.0 - *lam);
    for (int i = 0; i < *nt; i++) {
        double s = 0.0;
        for (int k = -*K; k <= *K; k++)
            s += line_density(theta[i] + 2.0 * M_PI * k, *mu, sigma, *beta,
                              *lam, 1.0 - *lam, c);
        out[i] = s;
    }
}

/* g_c(w) = int_0^inf cos(wu) exp(-u^beta) du and g_s (sine analogue), by composite Simpson on
 * [0, umax] with umax = 37^(1/beta) (exp(-37) below double resolution against the integral scale)
 * and a panel count resolving the oscillation. Both integrals are absolutely convergent for every
 * beta > 0. */
static void gcs(double w, double beta, double *gc, double *gs)
{
    double aw = fabs(w), sc = 0.0, ss = 0.0;
    /* [0, 1] under u = s^4, which removes the derivative singularity of exp(-u^beta) at zero */
    long n1 = 2000 + (long)(8.0 * aw / M_PI) * 2;
    if (n1 > 400000) n1 = 400000;
    double h1 = 1.0 / (double)n1;
    for (long j = 0; j <= n1; j++) {
        double s = h1 * (double)j;
        double u = s * s * s * s;
        double f = 4.0 * s * s * s * exp(-pow(u, beta));
        double wgt = (j == 0 || j == n1) ? 1.0 : ((j % 2 == 1) ? 4.0 : 2.0);
        sc += wgt * f * cos(aw * u);
        ss += wgt * f * sin(aw * u);
    }
    sc *= h1 / 3.0; ss *= h1 / 3.0;
    /* [1, umax] by plain composite Simpson; the integrand is smooth there */
    double umax = pow(37.0, 1.0 / beta);
    if (umax > 1.0) {
        double sc2 = 0.0, ss2 = 0.0;
        long n2 = 4000 + (long)(8.0 * aw * (umax - 1.0) / M_PI) * 2;
        if (n2 > 2000000) n2 = 2000000;
        double h2 = (umax - 1.0) / (double)n2;
        for (long j = 0; j <= n2; j++) {
            double u = 1.0 + h2 * (double)j;
            double f = exp(-pow(u, beta));
            double wgt = (j == 0 || j == n2) ? 1.0 : ((j % 2 == 1) ? 4.0 : 2.0);
            sc2 += wgt * f * cos(aw * u);
            ss2 += wgt * f * sin(aw * u);
        }
        sc += sc2 * h2 / 3.0; ss += ss2 * h2 / 3.0;
    }
    *gc = sc;
    *gs = (w < 0.0 ? -ss : ss);
}

void wsgn_cf(const double *t, const int *nt, const double *mu, const double *sigma,
             const double *beta, const double *lam1, const double *lam2,
             double *re, double *im)
{
    double cn = 1.0 / ((1.0 / *lam1 + 1.0 / *lam2) * tgamma((*beta + 1.0) / *beta));
    for (int i = 0; i < *nt; i++) {
        double g1c, g1s, g2c, g2s;
        gcs(*sigma * t[i] / *lam1, *beta, &g1c, &g1s);
        gcs(*sigma * t[i] / *lam2, *beta, &g2c, &g2s);
        double a = cn * (g1c / *lam1 + g2c / *lam2);          /* real part before rotation */
        double b = cn * (g1s / *lam1 - g2s / *lam2);          /* imaginary part before rotation */
        double cm = cos(*mu * t[i]), sm = sin(*mu * t[i]);
        re[i] = cm * a - sm * b;
        im[i] = sm * a + cm * b;
    }
}

void wsgn_trigmom(const int *P, const double *mu, const double *beta, const double *lam,
                  double *A, double *B)
{
    double sigma = sqrt(*lam * (1.0 - *lam));
    double l2 = 1.0 - *lam;
    for (int p = 1; p <= *P; p++) {
        double tp = (double)p;
        int one = 1;
        wsgn_cf(&tp, &one, mu, &sigma, beta, lam, &l2, &A[p - 1], &B[p - 1]);
    }
}

void wsgn_nll(const double *theta, const int *nt, const double *mu, const double *beta,
              const double *lam, const int *K, double *out)
{
    double sigma = sqrt(*lam * (1.0 - *lam));
    double c = norm_const(sigma, *beta, *lam, 1.0 - *lam);
    double s = 0.0;
    for (int i = 0; i < *nt; i++) {
        double f = 0.0;
        for (int k = -*K; k <= *K; k++)
            f += line_density(theta[i] + 2.0 * M_PI * k, *mu, sigma, *beta,
                              *lam, 1.0 - *lam, c);
        s -= log(f);
    }
    *out = s;
}

/* Regularised incomplete gamma functions P(a, x) and Q(a, x) = 1 - P(a, x), a > 0, x >= 0: the
 * power series where it converges fast (x < a + 1) and the Legendre continued fraction (modified
 * Lentz) otherwise, each returning the smaller tail directly so that neither is formed as a
 * difference of numbers near one. */
static void incgamma(double a, double x, double *P, double *Q)
{
    if (isnan(x)) { *P = NAN; *Q = NAN; return; }
    if (x <= 0.0) { *P = 0.0; *Q = 1.0; return; }
    if (isinf(x)) { *P = 1.0; *Q = 0.0; return; }
    double lpre = a * log(x) - x - lgamma(a);
    if (x < a + 1.0) {
        double ap = a, del = 1.0 / a, sum = del;
        for (int n = 0; n < 100000; n++) {
            ap += 1.0; del *= x / ap; sum += del;
            if (fabs(del) < fabs(sum) * 1e-16) break;
        }
        *P = sum * exp(lpre);
        if (*P > 1.0) *P = 1.0;
        *Q = 1.0 - *P;
    } else {
        const double tiny = 1e-300;
        double b = x + 1.0 - a, c = 1.0 / tiny, d = 1.0 / b, h = d;
        for (int i = 1; i < 100000; i++) {
            double an = -i * (i - a);
            b += 2.0;
            d = an * d + b; if (fabs(d) < tiny) d = tiny;
            c = b + an / c; if (fabs(c) < tiny) c = tiny;
            d = 1.0 / d;
            double del = d * c;
            h *= del;
            if (fabs(del - 1.0) < 1e-16) break;
        }
        *Q = exp(lpre) * h;
        if (*Q > 1.0) *Q = 1.0;
        *P = 1.0 - *Q;
    }
}

/* P(a < X <= b) for the line law, a <= b, from the two-piece incomplete-gamma form:
 * left of mu the mass beyond distance d is wL Q(1/beta, (lam2 d / sigma)^beta) with
 * wL = lam1 / (lam1 + lam2); right of mu it is wR Q(1/beta, (lam1 d / sigma)^beta) with
 * wR = lam2 / (lam1 + lam2). */
static double line_mass(double a, double b, double mu, double sigma, double beta,
                        double lam1, double lam2)
{
    double wL = lam1 / (lam1 + lam2), wR = lam2 / (lam1 + lam2), sh = 1.0 / beta;
    double Pa, Qa, Pb, Qb;
    if (b <= mu) {
        incgamma(sh, pow(lam2 * (mu - a) / sigma, beta), &Pa, &Qa);
        incgamma(sh, pow(lam2 * (mu - b) / sigma, beta), &Pb, &Qb);
        return wL * (Pa - Pb);
    }
    if (a >= mu) {
        incgamma(sh, pow(lam1 * (a - mu) / sigma, beta), &Pa, &Qa);
        incgamma(sh, pow(lam1 * (b - mu) / sigma, beta), &Pb, &Qb);
        return wR * (Pb - Pa);
    }
    incgamma(sh, pow(lam2 * (mu - a) / sigma, beta), &Pa, &Qa);
    incgamma(sh, pow(lam1 * (b - mu) / sigma, beta), &Pb, &Qb);
    return wL * Pa + wR * Pb;
}

void wsgn_pdist(const double *x, const int *nx, const double *mu, const double *sigma,
                const double *beta, const double *lam1, const double *lam2, double *out)
{
    double wL = *lam1 / (*lam1 + *lam2), wR = *lam2 / (*lam1 + *lam2), sh = 1.0 / *beta;
    for (int i = 0; i < *nx; i++) {
        double P, Q, d = x[i] - *mu;
        if (d < 0.0) {
            incgamma(sh, pow(*lam2 * (-d) / *sigma, *beta), &P, &Q);
            out[i] = wL * Q;
        } else {
            incgamma(sh, pow(*lam1 * d / *sigma, *beta), &P, &Q);
            out[i] = wL + wR * P;
        }
    }
}

/* wrapped distribution function on [0, 2*pi): the mass of the arcs [2 pi k, theta + 2 pi k],
 * |k| <= K, so that it is exactly the integral of wsgn_wdensity with the same K; at theta = 2 pi it
 * returns the mass the truncated wrapping sum captures, which is how the truncation is checked */
static double wcdf1(double theta, double mu, double beta, double lam, int K)
{
    double sigma = sqrt(lam * (1.0 - lam)), s = 0.0;
    for (int k = -K; k <= K; k++)
        s += line_mass(2.0 * M_PI * k, theta + 2.0 * M_PI * k, mu, sigma, beta, lam, 1.0 - lam);
    return s;
}

void wsgn_wcdf(const double *theta, const int *nt, const double *mu, const double *beta,
               const double *lam, const int *K, double *out)
{
    for (int i = 0; i < *nt; i++) {
        double t = theta[i];
        if (t <= 0.0) out[i] = 0.0;
        else out[i] = wcdf1(t < 2.0 * M_PI ? t : 2.0 * M_PI, *mu, *beta, *lam, *K);
    }
}

/* wrapped quantile function on [0, 2*pi), by bisection on the monotone distribution function;
 * 100 halvings take the bracket below the spacing of doubles near 2 pi */
void wsgn_wquant(const double *p, const int *np, const double *mu, const double *beta,
                 const double *lam, const int *K, double *out)
{
    for (int i = 0; i < *np; i++) {
        double lo = 0.0, hi = 2.0 * M_PI;
        if (p[i] <= 0.0) { out[i] = 0.0; continue; }
        if (p[i] >= 1.0) { out[i] = 2.0 * M_PI; continue; }
        for (int it = 0; it < 100; it++) {
            double mid = 0.5 * (lo + hi);
            if (mid <= lo || mid >= hi) break;
            if (wcdf1(mid, *mu, *beta, *lam, *K) < p[i]) lo = mid; else hi = mid;
        }
        out[i] = 0.5 * (lo + hi);
    }
}
