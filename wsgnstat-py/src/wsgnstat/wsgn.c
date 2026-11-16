/* wsgn.c : shared pure-C back-end for the two-piece exponential power distribution on the line and
 * its wrapped version on the circle (the package keeps the authors' earlier name, wrapped skew
 * generalised normal). See wsgn.h for the interface. Pure libm, so the same source compiles under R
 * (via .C) and Python (via ctypes) with byte-identical results. No random numbers are drawn here
 * except a fixed-seed splitmix64 stream for the starting values of grouped fits, which is identical
 * on every platform; samples for the parametric bootstrap are drawn by the front ends.
 *
 * Parameterisation of the wrapped family: location mu, scale sigma > 0, skewness 0 < lam < 1 and
 * shape beta > 0, with right rate a = lam / sigma and left rate b = (1 - lam) / sigma, so that the
 * line density is c exp{-(a (x - mu))^beta} for x >= mu and c exp{-(b (mu - x))^beta} for x < mu,
 * with c = lam (1 - lam) / {sigma Gamma(1 + 1/beta)} and P(X < mu) = lam. */
/* declare the reentrant lgamma_r (see wsgn_lgamma below) before <math.h> is first included */
#if defined(__APPLE__) && !defined(_REENTRANT)
#define _REENTRANT 1
#elif defined(__linux__) && !defined(_DEFAULT_SOURCE)
#define _DEFAULT_SOURCE 1
#endif
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "wsgn.h"

/* lgamma() writes the global signgam, so concurrent fits would race on it (harmlessly, but still a
 * race). The reentrant lgamma_r returns the same value where it is declared; elsewhere (Windows,
 * strict ISO C) the plain function is used. */
#if (defined(__APPLE__) && defined(_REENTRANT)) || (defined(__GLIBC__) && defined(_DEFAULT_SOURCE))
static double wsgn_lgamma(double x) { int sg; return lgamma_r(x, &sg); }
#else
static double wsgn_lgamma(double x) { return lgamma(x); }
#endif

#ifdef _OPENMP
#include <omp.h>
#endif
/* the parallel trials inside a fit use OpenMP tasks (taskloop, OpenMP 4.5); an older OpenMP keeps
 * only the parallel bootstrap */
#if defined(_OPENMP) && _OPENMP >= 201511
#define WSGN_TASKS 1
#define WSGN_PRAGMA(x) _Pragma(#x)
/* A loop of independent trials, run as OpenMP tasks when `par` is nonzero (which it is only inside
 * the parallel region opened by wsgn_fit or wsgn_profile) and otherwise in order on the calling
 * thread. Every iteration writes only to its own slot of the result arrays, and the caller combines
 * the slots in index order afterwards, so the result does not depend on the number of threads. */
#define WSGN_TASKLOOP(par) WSGN_PRAGMA(omp taskloop default(shared) grainsize(1) untied if(par))
#else
#define WSGN_TASKLOOP(par)
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define TWO_PI (2.0 * M_PI)
#define BIGOBJ 1e10
#define NOUT 12

/* 20-point Gauss-Legendre rule on [-1, 1] */
static const double GLX[20] = {
    0.99312859918509492479, 0.96397192727791379127, 0.91223442825132590587,
    0.83911697182221882339, 0.74633190646015079261, 0.63605368072651502545,
    0.510867001950827098, 0.37370608871541956067, 0.22778585114164507808,
    0.076526521133497333755, -0.076526521133497333755, -0.22778585114164507808,
    -0.37370608871541956067, -0.510867001950827098, -0.63605368072651502545,
    -0.74633190646015079261, -0.83911697182221882339, -0.91223442825132590587,
    -0.96397192727791379127, -0.99312859918509492479};
static const double GLW[20] = {
    0.017614007139152118312, 0.040601429800386941331, 0.06267204833410906357,
    0.083276741576704748725, 0.10193011981724043504, 0.11819453196151841731,
    0.1316886384491766269, 0.14209610931838205133, 0.14917298647260374679,
    0.1527533871307258507, 0.1527533871307258507, 0.14917298647260374679,
    0.14209610931838205133, 0.1316886384491766269, 0.11819453196151841731,
    0.10193011981724043504, 0.083276741576704748725, 0.06267204833410906357,
    0.040601429800386941331, 0.017614007139152118312};

static double logaddexp(double a, double b)
{
    if (a == -INFINITY) return b;
    if (b == -INFINITY) return a;
    return (a > b) ? a + log1p(exp(b - a)) : b + log1p(exp(a - b));
}

/* angle reduced to [0, 2 pi) */
static double mod2pi(double x)
{
    double t = fmod(x, TWO_PI);
    if (t < 0.0) t += TWO_PI;
    if (t >= TWO_PI) t -= TWO_PI;
    return t;
}

/* ------------------------------------------------------------------ incomplete gamma */

/* Regularised incomplete gamma functions P(a, x), Q(a, x) = 1 - P(a, x) and log Q(a, x), a > 0,
 * x >= 0: the power series where it converges fast (x < a + 1) and the Legendre continued fraction
 * (modified Lentz) otherwise, each returning the smaller tail directly, so that neither is formed
 * as a difference of numbers near one, and log Q without underflow far in the tail. */
/* core with both x and lx = log x, so that x may underflow to 0 while a log x stays exact (the
 * lower tail of a shape far above one next to the location) */
static void incgamma_core(double a, double x, double lx, double *P, double *Q, double *lQ)
{
    if (isnan(x) || isnan(lx)) { *P = NAN; *Q = NAN; *lQ = NAN; return; }
    if (lx == -INFINITY || x < 0.0) { *P = 0.0; *Q = 1.0; *lQ = 0.0; return; }
    if (isinf(x)) { *P = 1.0; *Q = 0.0; *lQ = -INFINITY; return; }
    double lpre = a * lx - x - wsgn_lgamma(a);
    if (x < a + 1.0) {
        double ap = a, del = 1.0 / a, sum = del;
        for (int n = 0; n < 100000; n++) {
            ap += 1.0; del *= x / ap; sum += del;
            if (fabs(del) < fabs(sum) * 1e-17) break;
        }
        *P = sum * exp(lpre);
        if (*P > 1.0) *P = 1.0;
        *Q = 1.0 - *P;
        *lQ = log1p(-*P);
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
            if (fabs(del - 1.0) < 4.5e-16) break;   /* 2 ulp: 1e-16 can never be met */
        }
        *lQ = lpre + log(h);
        *Q = exp(*lQ);
        if (*Q > 1.0) { *Q = 1.0; *lQ = 0.0; }
        *P = 1.0 - *Q;
    }
}

static void incgamma3(double a, double x, double *P, double *Q, double *lQ)
{
    incgamma_core(a, x, (x > 0.0) ? log(x) : (x == 0.0 ? -INFINITY : NAN), P, Q, lQ);
}

/* P(a, x) and Q(a, x) from lx = log x */
static void incgamma_lx(double a, double lx, double *P, double *Q)
{
    double lQ;
    incgamma_core(a, exp(lx), lx, P, Q, &lQ);
}

static double log_incgamma_Q(double a, double x)
{
    double P, Q, lQ;
    incgamma3(a, x, &P, &Q, &lQ);
    return lQ;
}

/* ------------------------------------------------------------------ the line distribution */

/* normalising constant of the general two-piece line density, right rate lam1 / sigma and left
 * rate lam2 / sigma */
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

void wsgn_pdist(const double *x, const int *nx, const double *mu, const double *sigma,
                const double *beta, const double *lam1, const double *lam2, double *out)
{
    double wL = *lam1 / (*lam1 + *lam2), wR = *lam2 / (*lam1 + *lam2), sh = 1.0 / *beta;
    for (int i = 0; i < *nx; i++) {
        double P, Q, d = x[i] - *mu;
        if (isnan(d)) { out[i] = NAN; continue; }
        if (d < 0.0) {
            incgamma_lx(sh, *beta * log(*lam2 * (-d) / *sigma), &P, &Q);
            out[i] = wL * Q;
        } else {
            incgamma_lx(sh, (d > 0.0) ? *beta * log(*lam1 * d / *sigma) : -INFINITY, &P, &Q);
            out[i] = wL + wR * P;
        }
    }
}

/* ------------------------------------------------------------------ characteristic function */

/* Gn(w) = G(w) / Gamma(1 + 1/beta), G(w) = int_0^inf exp(i w u - u^beta) du, w > 0, by rotating
 * the path of integration to the ray u = r exp(i phi). Along the arc joining the real axis to the
 * ray, |exp(i w u)| <= 1 and exp(-u^beta) decays as long as beta * phi < pi / 2, so the integral is
 * unchanged for every 0 <= phi < min(pi, pi / (2 beta)); phi is taken in the middle of that range.
 * On the ray both factors decay and the oscillation is slow, and the substitution r = exp(x) turns
 * the half-line into the whole real line with double-exponential decay at the right end and
 * exponential decay at the left. The integrand is analytic in a strip about the real x-axis, so the
 * trapezoid rule converges exponentially in 1/h; h is halved until two successive sums agree to
 * tol relative to the L1 norm of the integrand. This removes the long range [0, 37^(1/beta)] that
 * defeats direct quadrature for small shapes. */
static void G_quad(double w, double beta, double tol, double *gr, double *gi)
{
    if (w == 0.0) { *gr = 1.0; *gi = 0.0; return; }
    int neg = (w < 0.0);
    double aw = fabs(w);
    double phimax = (beta <= 0.5) ? M_PI : M_PI / (2.0 * beta);
    double phi = 0.5 * phimax;
    double sp = sin(phi), cp = cos(phi), cb = cos(beta * phi), sb = sin(beta * phi);
    double lgam = wsgn_lgamma(1.0 + 1.0 / beta);
    double lS0 = -log(aw) - lgam;
    if (lS0 > 0.0) lS0 = 0.0;
    double target = log(1e-3 * tol) + lS0;
    /* left end: |integrand| <= exp(x - lgam) */
    double xlo = target + lgam;
    /* right end: step out until the log-modulus is below target and falling at rate at least 1/2 */
    double xhi = xlo;
    for (int it = 0; it < 100000; it++) {
        double r = exp(xhi), rb = exp(beta * xhi);
        double lm = xhi - aw * r * sp - rb * cb - lgam;
        double slope = 1.0 - aw * r * sp - beta * rb * cb;
        if (lm < target && slope < -0.5) break;
        xhi += 0.25;
    }
    int N = 64;
    double h = (xhi - xlo) / (double)N;
    double sr = 0.0, si = 0.0, sa = 0.0;
    for (int j = 0; j <= N; j++) {
        double x = xlo + h * (double)j;
        double r = exp(x), rb = exp(beta * x);
        double ere = -aw * r * sp - rb * cb;
        double eim = aw * r * cp - rb * sb + phi;
        double m = exp(x + ere - lgam);
        double wt = (j == 0 || j == N) ? 0.5 : 1.0;
        sr += wt * m * cos(eim); si += wt * m * sin(eim); sa += wt * m;
    }
    double Tr = h * sr, Ti = h * si;
    for (int lev = 0; lev < 16; lev++) {
        h *= 0.5;
        for (int j = 1; j < 2 * N; j += 2) {
            double x = xlo + h * (double)j;
            double r = exp(x), rb = exp(beta * x);
            double ere = -aw * r * sp - rb * cb;
            double eim = aw * r * cp - rb * sb + phi;
            double m = exp(x + ere - lgam);
            sr += m * cos(eim); si += m * sin(eim); sa += m;
        }
        N *= 2;
        double Nr = h * sr, Ni = h * si, L1 = h * sa;
        double dr = Nr - Tr, di = Ni - Ti;
        Tr = Nr; Ti = Ni;
        if (lev >= 2 && sqrt(dr * dr + di * di) <= tol * L1) break;
    }
    *gr = Tr;
    *gi = neg ? -Ti : Ti;
}

/* The same Gn(w) from its series (the residue series of its Mellin-Barnes form): the descending series
 * sum_k (-1)^k Gamma(k beta + 1) exp(i pi (k beta + 1) / 2) w^-(k beta + 1) / k!, convergent for
 * every w > 0 when beta < 1, and the ascending series (1/beta) sum_k Gamma((k + 1)/beta) (i w)^k / k!,
 * convergent for every w when beta > 1; at beta = 1, G(w) = 1 / (1 - i w). Both series cancel badly
 * in parts of the (w, beta) plane. The rounding error is estimated, with a safety factor of about 4
 * over the unit roundoff, as 1e-15 sum_k |t_k| (1 + |lt_k| + |angle_k|), lt_k being the log-gamma combination whose exponential is the term, and NaN is
 * returned when that exceeds tol |sum|, or when the sum breaks |Gn| <= 1. Kept as a check on
 * G_quad. */
static void G_series(double w, double beta, double tol, double *gr, double *gi)
{
    if (w == 0.0) { *gr = 1.0; *gi = 0.0; return; }
    int neg = (w < 0.0);
    double aw = fabs(w);
    double sr = 0.0, si = 0.0, tmax = 0.0, rerr = 0.0;
    if (beta == 1.0) {
        sr = 1.0 / (1.0 + aw * aw); si = aw / (1.0 + aw * aw);
        *gr = sr; *gi = neg ? -si : si;
        return;
    }
    double lgam = wsgn_lgamma(1.0 + 1.0 / beta);
    int small = 0;
    int conv = 0;
    for (int k = 0; k < 200000; k++) {
        double lt, ang;
        if (beta < 1.0) {
            lt = wsgn_lgamma(k * beta + 1.0) - (k * beta + 1.0) * log(aw) - wsgn_lgamma(k + 1.0) - lgam;
            ang = M_PI * (k * beta + 1.0) / 2.0 + ((k % 2) ? M_PI : 0.0);
        } else {
            lt = wsgn_lgamma((k + 1.0) / beta) - wsgn_lgamma(1.0 / beta) + k * log(aw) - wsgn_lgamma(k + 1.0);
            ang = M_PI * (double)(k % 4) / 2.0;
        }
        double t = exp(lt);
        if (t > tmax) tmax = t;
        sr += t * cos(ang); si += t * sin(ang);
        /* rounding error of this term: exp of a log-gamma combination of size |lt|, and an angle */
        rerr += t * (1.0 + fabs(lt) + fabs(ang));
        /* since |Gn| <= 1, a rounding error above tol already means refusal */
        if (!(1e-15 * rerr <= tol)) { tmax = INFINITY; break; }
        double s = sqrt(sr * sr + si * si);
        /* stop once past the largest term and the terms are negligible against both the largest
         * term and the tolerance on the current sum */
        if (k > 2 && t < tmax && t < 1e-20 * tmax && t < 1e-3 * tol * s) {
            if (++small >= 3) { conv = 1; break; }
        } else small = 0;
    }
    double s = sqrt(sr * sr + si * si);
    if (!conv) rerr = INFINITY;
    /* refuse when the series has not converged within 200000 terms, when the accumulated
     * rounding error exceeds tol relative to the sum, when anything is not finite, or when the sum
     * breaks the bound |Gn| <= 1 */
    if (!(isfinite(rerr) && isfinite(s) && s <= 1.0 + 1e-12 && 1e-15 * rerr <= tol * s)) {
        sr = NAN; si = NAN;
    }
    *gr = sr; *gi = neg ? -si : si;
}

void wsgn_cf(const double *t, const int *nt, const double *mu, const double *sigma,
             const double *beta, const double *lam1, const double *lam2, const double *tol0,
             const int *method, double *re, double *im)
{
    /* below 1e-14 the stopping test of the quadrature is beyond rounding; the fronts raise smaller
     * tolerances to this floor with a warning, and so does the C code */
    double tolv = (*tol0 >= 1e-14) ? *tol0 : 1e-14;
    const double *tol = &tolv;
    /* phi(t) = exp(i t mu) [w1 Gn(t / a1) + w2 conj Gn(t / a2)], a1 = lam1 / sigma, a2 =
     * lam2 / sigma, with w1 = P(X > mu) = lam2 / (lam1 + lam2) and w2 = P(X < mu) */
    double w1 = *lam2 / (*lam1 + *lam2), w2 = *lam1 / (*lam1 + *lam2);
    for (int i = 0; i < *nt; i++) {
        double g1r, g1i, g2r, g2i;
        if (*method == 1) {
            G_series(*sigma * t[i] / *lam1, *beta, *tol, &g1r, &g1i);
            G_series(*sigma * t[i] / *lam2, *beta, *tol, &g2r, &g2i);
        } else {
            G_quad(*sigma * t[i] / *lam1, *beta, *tol, &g1r, &g1i);
            G_quad(*sigma * t[i] / *lam2, *beta, *tol, &g2r, &g2i);
        }
        double a = w1 * g1r + w2 * g2r;
        double b = w1 * g1i - w2 * g2i;
        double cm = cos(*mu * t[i]), sm = sin(*mu * t[i]);
        re[i] = cm * a - sm * b;
        im[i] = sm * a + cm * b;
    }
}

void wsgn_trigmom(const int *P, const double *mu, const double *sigma, const double *beta,
                  const double *lam, const double *tol, double *A, double *B)
{
    double l2 = 1.0 - *lam;
    int one = 1, meth = 0;
    for (int p = 1; p <= *P; p++) {
        double tp = (double)p;
        wsgn_cf(&tp, &one, mu, sigma, beta, lam, &l2, tol, &meth, &A[p - 1], &B[p - 1]);
    }
}

void wsgn_decay(const double *sigma, const double *beta, const double *lam, double *out)
{
    /* from the first two terms of the descending series: |phi_p| ~ C p^-(1 + beta), with
     * C = c Gamma(beta + 1) |a^beta exp(i pi (beta + 1)/2) + b^beta exp(-i pi (beta + 1)/2)| */
    double be = *beta, la = log(*lam / *sigma), lb = log((1.0 - *lam) / *sigma);
    double lc = log(*lam) + log1p(-*lam) - log(*sigma) - wsgn_lgamma(1.0 + 1.0 / be);
    double A = exp(be * la) + exp(be * lb), D = exp(be * la) - exp(be * lb);
    double s = sin(M_PI * be / 2.0), co = cos(M_PI * be / 2.0);
    if (be == 2.0 * floor(be / 2.0)) s = 0.0;   /* sin(pi beta / 2) is exactly zero at even integers */
    out[0] = exp(lc + wsgn_lgamma(be + 1.0)) * sqrt(A * A * s * s + D * D * co * co);
    out[1] = 1.0 + be;
}

/* ------------------------------------------------------------------ wrapped density */

/* B_2j / (2j)! for j = 1, ..., 7 */
static const double BERN[7] = {1.0 / 12.0, -1.0 / 720.0, 1.0 / 30240.0, -1.0 / 1209600.0,
                               1.0 / 47900160.0, -691.0 / 1307674368000.0, 1.0 / 74724249600.0};

/* Euler-Maclaurin pieces for f(x) = exp{-(rate x)^beta} at x with w = (rate x)^beta:
 * sum_{k >= 0} f(x + h k) = (1/h) int_x^inf f + f(x) corr, with
 * corr = 1/2 - sum_{j=1}^{6} B_2j / (2j)! h^(2j-1) f^(2j-1)(x) / f(x), and err the size, relative to
 * f(x), of the first omitted term (j = 7). With g_n = d^n log f / dx^n
 * = -w beta (beta - 1) ... (beta - n + 1) / x^n, the ratios r_n = f^(n) / f follow from
 * r_(n+1) = sum_k C(n, k) g_(k+1) r_(n-k). Returns 0 when the formula is not usable at x (the
 * summand falls by more than a factor exp(1/2) per step, or corr is not positive). */
static int em_terms(double x, double w, double beta, double h, double *corr, double *err)
{
    double g[15], r[15];
    if (!(h * beta * w / x <= 0.5)) return 0;
    g[1] = -beta * w / x;
    for (int n = 2; n <= 14; n++) g[n] = g[n - 1] * (beta - (double)(n - 1)) / x;
    r[0] = 1.0;
    for (int n = 0; n < 13; n++) {
        double s = 0.0, binom = 1.0;
        for (int k = 0; k <= n; k++) {
            s += binom * g[k + 1] * r[n - k];
            binom = binom * (double)(n - k) / (double)(k + 1);
        }
        r[n + 1] = s;
    }
    double c = 0.5, hp = h;
    for (int j = 1; j <= 6; j++) {
        c -= BERN[j - 1] * hp * r[2 * j - 1];
        hp *= h * h;
    }
    *corr = c;
    *err = fabs(BERN[6] * hp * r[13]);
    return isfinite(c) && isfinite(*err) && c > 0.0;
}

/* constants of one side of the wrapped density: log rate, shape, 1/shape, wsgn_lgamma(1/shape), log shape */
typedef struct { double lrate, beta, alpha, lga, lbe; } sidec;

static sidec side_const(double lrate, double beta)
{
    sidec c;
    c.lrate = lrate; c.beta = beta; c.alpha = 1.0 / beta;
    c.lga = wsgn_lgamma(c.alpha); c.lbe = log(beta);
    return c;
}

/* Try to finish a side sum from the image at xk (w = (rate xk)^beta) on: returns 1 and the log of
 * the whole sum (in units of exp(-w0)) in *res if the rest is negligible or the Euler-Maclaurin tail
 * is accurate there (first omitted term below 1e-16 of the sum). Everything is in log space, since
 * the tail integral carries Gamma(1/beta), which overflows for beta below about 1/171. */
static int side_finish(double xk, double w, double w0, double ls, const sidec *sc, double *res)
{
    const double h = TWO_PI, lh = log(TWO_PI);
    double lfk = -(w - w0);
    double lint = sc->lga + log_incgamma_Q(sc->alpha, w) - lh - sc->lbe - sc->lrate + w0;
    if (logaddexp(lfk, lint) < ls - 39.2) { *res = ls; return 1; }
    double corr, err;
    if (em_terms(xk, w, sc->beta, h, &corr, &err)) {
        double ltot = logaddexp(ls, lint);
        if (log(err) + lfk <= ltot - 36.9) {
            *res = logaddexp(ltot, lfk + log(corr));
            return 1;
        }
    }
    return 0;
}

/* log of sum_{k >= 0} exp{-(rate (d0 + 2 pi k))^beta}, d0 >= 0 (the factor c is added by the
 * caller). Images are summed directly; after 8, 16 and 32 of them, after K of them, and after each
 * further image, the sum is finished as soon as the rest is negligible or the Euler-Maclaurin tail
 * (exact incomplete-gamma integral, Bernoulli terms through B_12) is accurate. The result does not
 * depend on K beyond rounding; K = 50 by default bounds the direct sum when the tail is not yet
 * usable. */
static double side_logsum(double d0, const sidec *sc, int K)
{
    const double h = TWO_PI;
    double beta = sc->beta, lrate = sc->lrate, res;
    if (isnan(d0) || isnan(lrate) || isnan(beta)) return NAN;
    double w0 = (d0 > 0.0) ? exp(beta * (lrate + log(d0))) : 0.0;
    if (isinf(w0)) return -INFINITY;   /* (rate d0)^beta overflows: this side contributes nothing */
    double s = 0.0, next_check = 30.0;
    int k;
    for (k = 0; k < K + 1000000; k++) {
        double x = d0 + h * (double)k;
        double w = (x > 0.0) ? exp(beta * (lrate + log(x))) : 0.0;
        if (isinf(w)) return -w0 + log(s);   /* every further image underflows */
        if (k > 0 && (k == 8 || k == 16 || k == 32 || k >= K) &&
            side_finish(x, w, w0, log(s), sc, &res))
            return -w0 + res;
        double dw = w - w0;
        s += exp(-dw);
        if (dw > next_check) {
            double lrem = sc->lga + log_incgamma_Q(sc->alpha, w) - log(h) - sc->lbe - lrate;
            if (lrem < -w0 + log(s) - 36.0) return -w0 + log(s);
            next_check = dw + 3.0;
        }
    }
    return -w0 + log(s);
}

/* the wrapped density's constants for one parameter value */
typedef struct { double mu, lc; sidec R, L; } wconst;

static wconst wconst_rates(double mu, double la, double lb, double beta)
{
    wconst w;
    w.mu = mu;
    w.lc = log(beta) - wsgn_lgamma(1.0 / beta) - logaddexp(-la, -lb);
    w.R = side_const(la, beta);
    w.L = side_const(lb, beta);
    return w;
}

static double logdens_w(double theta, const wconst *w, int K)
{
    /* the two distances from mu sum to 2 pi. The signed difference is reduced by remainder(),
     * which is exact, and the side is read from its sign; the distance on the side the point lies
     * is then the exact difference, and the other is 2 pi minus it. (Reducing mu - theta separately
     * with a modulus would round a tiny negative distance up to 2 pi and put the point at
     * distance 0 on the wrong side of mu.) */
    double d = remainder(theta - w->mu, TWO_PI), t, tl;
    if (d >= 0.0) { t = d; tl = TWO_PI - d; } else { tl = -d; t = TWO_PI + d; }
    return w->lc + logaddexp(side_logsum(t, &w->R, K), side_logsum(tl, &w->L, K));
}

/* log wrapped density at theta, rates on the log scale */
static double logdens_rates(double theta, double mu, double la, double lb, double beta, int K)
{
    wconst w = wconst_rates(mu, la, lb, beta);
    return logdens_w(theta, &w, K);
}

void wsgn_wdensity(const double *theta, const int *nt, const double *mu, const double *sigma,
                   const double *beta, const double *lam, const int *K, const int *give_log,
                   double *out)
{
    wconst w = wconst_rates(*mu, log(*lam / *sigma), log((1.0 - *lam) / *sigma), *beta);
    for (int i = 0; i < *nt; i++) {
        double v = logdens_w(theta[i], &w, *K);
        out[i] = *give_log ? v : exp(v);
    }
}

void wsgn_nll(const double *theta, const int *nt, const double *mu, const double *sigma,
              const double *beta, const double *lam, const int *K, double *out)
{
    wconst w = wconst_rates(*mu, log(*lam / *sigma), log((1.0 - *lam) / *sigma), *beta);
    double s = 0.0;
    for (int i = 0; i < *nt; i++)
        s -= logdens_w(theta[i], &w, *K);
    *out = s;
}

/* ------------------------------------------------------------------ arc probabilities */

/* sum_{k >= K} f(x + 2 pi k) for the side-distance density f(y) = (a beta / Gamma(1/beta))
 * exp{-(a y)^beta}: images are added directly until the Euler-Maclaurin tail is accurate (its
 * first omitted term below 1e-16 of ref, the density scale of the images already summed, plus the
 * tail) or the rest is negligible */
static double side_tail_density(double x, double lrate, double beta, int K, double ref)
{
    const double h = TWO_PI;
    double alpha = 1.0 / beta;
    double lf0 = lrate + log(beta) - wsgn_lgamma(alpha);
    double s = 0.0;
    for (int k = K; k < K + 1000000; k++) {
        double y = x + h * (double)k;
        double w = exp(beta * (lrate + log(y)));
        double fk = exp(lf0 - w);
        double corr, err;
        if (em_terms(y, w, beta, h, &corr, &err)) {
            double tint = exp(log_incgamma_Q(alpha, w) - log(h));
            if (err * fk <= 1e-16 * (ref + s + tint)) return s + tint + fk * corr;
        }
        s += fk;
        if (fk <= 1e-18 * (s + ref) || fk == 0.0) return s;
    }
    return s;
}

/* P(D in [x, x + wd]) for one image, the side distance D having density
 * f(y) = (a beta / Gamma(1/beta)) exp{-(a y)^beta}: from the incomplete gamma function when the
 * interval reaches back to near the origin or f varies strongly across it, otherwise (a short
 * interval on which log f varies by less than 1) by the 20-point Gauss-Legendre rule, so that a
 * short interval keeps its relative accuracy */
static double image_mass(double x, double wd, double lrate, double beta)
{
    double alpha = 1.0 / beta;
    double x1 = x + wd;
    double steep = beta * exp(beta * (lrate + log(x1))) * wd / (x > 0.0 ? x : 1.0);
    if (x <= wd || steep >= 1.0) {
        double P0, Q0, P1, Q1;
        incgamma_lx(alpha, (x > 0.0) ? beta * (lrate + log(x)) : -INFINITY, &P0, &Q0);
        incgamma_lx(alpha, beta * (lrate + log(x1)), &P1, &Q1);
        return (Q1 < 0.5) ? (Q0 - Q1) : (P1 - P0);
    }
    double lf0 = lrate + log(beta) - wsgn_lgamma(alpha), s = 0.0, r = 0.5 * wd;
    for (int i = 0; i < 20; i++) {
        double y = x + r * (1.0 + GLX[i]);
        s += GLW[i] * exp(lf0 - exp(beta * (lrate + log(y))));
    }
    return r * s;
}

/* sum_{k >= 0} P(D in [y + 2 pi k, y + wd + 2 pi k]), y >= 0, 0 <= wd <= 2 pi. Wide arcs: exact
 * incomplete-gamma differences for k < K (each formed from whichever of P and Q keeps it accurate),
 * stopping once the remaining mass is negligible, then the Gauss-Legendre integral of the
 * Euler-Maclaurin tail density. Short arcs (wd < 0.01): the first image by image_mass and all the
 * others through the tail density integrated across the arc. The width is passed separately from
 * the start, so a short arc keeps its relative accuracy. */
static double side_arcmass(double y, double wd, double lrate, double beta, int K)
{
    const double h = TWO_PI;
    if (!(wd > 0.0)) return 0.0;
    if (y < 0.0) y = 0.0;
    double alpha = 1.0 / beta, tot = 0.0, s = 0.0, r = 0.5 * wd;
    int k0 = 0;
    if (wd < 0.01) {
        tot = image_mass(y, wd, lrate, beta);
        k0 = 1;
    } else {
        for (int k = 0; k < K; k++) {
            double x0 = y + h * (double)k, x1 = x0 + wd;
            double P0, Q0, P1, Q1;
            incgamma_lx(alpha, (x0 > 0.0) ? beta * (lrate + log(x0)) : -INFINITY, &P0, &Q0);
            incgamma_lx(alpha, beta * (lrate + log(x1)), &P1, &Q1);
            tot += (Q1 < 0.5) ? (Q0 - Q1) : (P1 - P0);
            if (Q1 <= 1e-17 * tot || Q1 < 1e-300) return tot;
        }
        k0 = K;
    }
    for (int i = 0; i < 20; i++)
        s += GLW[i] * side_tail_density(y + r * (1.0 + GLX[i]), lrate, beta, k0, tot / wd);
    return tot + r * s;
}

/* P(Theta in the arc of width wd measured anticlockwise from lo) */
static double arc_prob(double lo, double wd, double mu, double la, double lb, double beta, int K)
{
    if (!(wd > 0.0)) return 0.0;
    if (wd >= TWO_PI) return 1.0;
    double pR = 1.0 / (1.0 + exp(la - lb)), pL = 1.0 / (1.0 + exp(lb - la));
    /* the distribution function is continuous in lo, so a start within rounding of mu on either
     * side changes the arc probability by at most that distance times the density: mod2pi is safe */
    double t0 = mod2pi(lo - mu);
    double su[2], sv[2];
    int ns = 1;
    if (t0 + wd <= TWO_PI) { su[0] = t0; sv[0] = wd; }
    else { su[0] = t0; sv[0] = TWO_PI - t0; su[1] = 0.0; sv[1] = wd - sv[0]; ns = 2; }
    double tot = 0.0;
    for (int s = 0; s < ns; s++) {
        tot += pR * side_arcmass(su[s], sv[s], la, beta, K);
        tot += pL * side_arcmass((TWO_PI - su[s]) - sv[s], sv[s], lb, beta, K);
    }
    return tot;
}

void wsgn_wcdf(const double *theta, const int *nt, const double *mu, const double *sigma,
               const double *beta, const double *lam, const int *K, double *out)
{
    double la = log(*lam / *sigma), lb = log((1.0 - *lam) / *sigma);
    for (int i = 0; i < *nt; i++) {
        double t = theta[i];
        if (isnan(t)) out[i] = NAN;
        else if (t <= 0.0) out[i] = 0.0;
        else if (t >= TWO_PI) out[i] = 1.0;
        else out[i] = arc_prob(0.0, t, *mu, la, lb, *beta, *K);
    }
}

/* wrapped quantile function on [0, 2 pi], by bisection on the monotone distribution function:
 * geometric midpoints while the bracket spans more than a factor 4, so that small probabilities
 * get quantiles with relative accuracy, then ordinary halving to the spacing of doubles */
void wsgn_wquant(const double *p, const int *np, const double *mu, const double *sigma,
                 const double *beta, const double *lam, const int *K, double *out)
{
    double la = log(*lam / *sigma), lb = log((1.0 - *lam) / *sigma);
    const double qmin = 1e-300;
    double pmin = arc_prob(0.0, qmin, *mu, la, lb, *beta, *K);
    for (int i = 0; i < *np; i++) {
        if (isnan(p[i])) { out[i] = NAN; continue; }
        if (p[i] <= 0.0) { out[i] = 0.0; continue; }
        if (p[i] >= 1.0) { out[i] = TWO_PI; continue; }
        if (p[i] <= pmin) { out[i] = (pmin > 0.0) ? qmin * (p[i] / pmin) : 0.0; continue; }
        double lo = qmin, hi = TWO_PI;
        for (int it = 0; it < 300; it++) {
            double mid = (hi > 4.0 * lo) ? sqrt(lo) * sqrt(hi) : 0.5 * (lo + hi);
            if (mid <= lo || mid >= hi) break;
            if (arc_prob(0.0, mid, *mu, la, lb, *beta, *K) < p[i]) lo = mid; else hi = mid;
        }
        out[i] = 0.5 * (lo + hi);
    }
}

/* ------------------------------------------------------------------ maximum likelihood */

typedef struct {
    const double *th; int n;                 /* exact observations (continuous likelihood) */
    const double *cnt; int J; double hc, orig;  /* grouped data: counts on a regular grid of cells */
    int grouped, sym, K;
    int fixmu_on; double fixmu;
    int fixbeta_on; double fixbeta;
    double beta0, beta1, rlim;
    int par;                                 /* nonzero: independent trials may run as OpenMP tasks */
} fctx;

static int ctx_npar(const fctx *c) { return 4 - c->fixmu_on - c->sym - c->fixbeta_on; }

static void unpack(const double *x, const fctx *c, double *mu, double *la, double *lb, double *lbe)
{
    int i = 0;
    *mu = c->fixmu_on ? c->fixmu : x[i++];
    *la = x[i++];
    *lb = c->sym ? *la : x[i++];
    *lbe = c->fixbeta_on ? log(c->fixbeta) : x[i++];
}

static void pack(double *x, const fctx *c, double mu, double la, double lb, double lbe)
{
    int i = 0;
    if (!c->fixmu_on) x[i++] = mu;
    x[i++] = la;
    if (!c->sym) x[i++] = lb;
    if (!c->fixbeta_on) x[i++] = lbe;
}

/* the terms of one likelihood evaluation computed as tasks in fixed chunks, then summed on the
 * calling thread in index order, so the value is bit-identical to the serial sum */
#define WSGN_CHUNK 64
static double loglik_ctx(const fctx *c, double mu, double la, double lb, double be)
{
    double v = 0.0;
#ifdef WSGN_TASKS
    int m = c->grouped ? c->J : c->n;
    if (c->par && m >= 2 * WSGN_CHUNK) {
        double *t = (double *)malloc(sizeof(double) * (size_t)m);
        if (t) {
            int nch = (m + WSGN_CHUNK - 1) / WSGN_CHUNK;
            if (c->grouped) {
#pragma omp taskloop default(shared) grainsize(1)
                for (int k = 0; k < nch; k++) {
                    int hi = (k + 1) * WSGN_CHUNK < m ? (k + 1) * WSGN_CHUNK : m;
                    for (int j = k * WSGN_CHUNK; j < hi; j++) {
                        if (c->cnt[j] <= 0.0) { t[j] = 0.0; continue; }
                        double ce = c->orig + c->hc * (double)j;
                        t[j] = arc_prob(ce - 0.5 * c->hc, c->hc, mu, la, lb, be, c->K);
                    }
                }
                for (int j = 0; j < m; j++) {
                    if (c->cnt[j] <= 0.0) continue;
                    double p = t[j];
                    if (!(p > 0.0)) { free(t); return -INFINITY; }
                    v += c->cnt[j] * log(p);
                }
            } else {
                wconst w = wconst_rates(mu, la, lb, be);
#pragma omp taskloop default(shared) grainsize(1)
                for (int k = 0; k < nch; k++) {
                    int hi = (k + 1) * WSGN_CHUNK < m ? (k + 1) * WSGN_CHUNK : m;
                    for (int i = k * WSGN_CHUNK; i < hi; i++) t[i] = logdens_w(c->th[i], &w, c->K);
                }
                for (int i = 0; i < m; i++) v += t[i];
            }
            free(t);
            return v;
        }
    }
#endif
    if (c->grouped) {
        for (int j = 0; j < c->J; j++) {
            if (c->cnt[j] <= 0.0) continue;
            double ce = c->orig + c->hc * (double)j;
            double p = arc_prob(ce - 0.5 * c->hc, c->hc, mu, la, lb, be, c->K);
            if (!(p > 0.0)) return -INFINITY;
            v += c->cnt[j] * log(p);
        }
    } else {
        wconst w = wconst_rates(mu, la, lb, be);
        for (int i = 0; i < c->n; i++) v += logdens_w(c->th[i], &w, c->K);
    }
    return v;
}

static double objective(const double *x, const fctx *c)
{
    double mu, la, lb, lbe;
    unpack(x, c, &mu, &la, &lb, &lbe);
    double be = exp(lbe);
    if (!c->fixbeta_on && !(be >= c->beta0 && be <= c->beta1)) return BIGOBJ;
    if (!(fabs(la) <= c->rlim && fabs(lb) <= c->rlim)) return BIGOBJ;
    double v = -loglik_ctx(c, mu, la, lb, be);
    return isfinite(v) ? v : BIGOBJ;
}

/* Nelder-Mead simplex, following scipy.optimize.minimize(method = "Nelder-Mead") with its default
 * coefficients, initial simplex and stopping rule (xatol, fatol); d <= 4 */
static double nelder_mead(const fctx *c, double *x, int d, int maxfev, double xatol, double fatol,
                          int *nfev)
{
    const double rho = 1.0, chi = 2.0, psi = 0.5, sig = 0.5;
    double sim[5][4], fs[5], xbar[4], xr[4], xe[4], xc[4], tmp[4];
    int fe = 0, it = 0, j, k;
    for (j = 0; j < d; j++) sim[0][j] = x[j];
    for (k = 0; k < d; k++) {
        for (j = 0; j < d; j++) sim[k + 1][j] = x[j];
        sim[k + 1][k] = (x[k] != 0.0) ? 1.05 * x[k] : 0.00025;
    }
    for (k = 0; k <= d; k++) { fs[k] = objective(sim[k], c); fe++; }
    /* sort */
    for (k = 1; k <= d; k++) {
        double fv = fs[k]; memcpy(tmp, sim[k], sizeof tmp);
        int m = k - 1;
        while (m >= 0 && fs[m] > fv) { fs[m + 1] = fs[m]; memcpy(sim[m + 1], sim[m], sizeof tmp); m--; }
        fs[m + 1] = fv; memcpy(sim[m + 1], tmp, sizeof tmp);
    }
    while (fe < maxfev && it < maxfev) {
        double xd = 0.0, fd = 0.0, xr0 = 0.0;
        for (k = 1; k <= d; k++) {
            for (j = 0; j < d; j++) {
                double a = fabs(sim[k][j] - sim[0][j]), ar = a / (1.0 + fabs(sim[0][j]));
                if (a > xd) xd = a;
                if (ar > xr0) xr0 = ar;
            }
            double b = fabs(fs[0] - fs[k]); if (b > fd) fd = b;
        }
        if (xd <= xatol && fd <= fatol) break;
        /* a simplex collapsed to a few units in the last place in every coordinate cannot improve; the function
         * tolerance may never be met next to a cusp, where the likelihood changes by orders of
         * magnitude between neighbouring doubles of the location */
        if (xr0 <= 2e-15) break;
        for (j = 0; j < d; j++) {
            xbar[j] = 0.0;
            for (k = 0; k < d; k++) xbar[j] += sim[k][j];
            xbar[j] /= (double)d;
        }
        for (j = 0; j < d; j++) xr[j] = (1.0 + rho) * xbar[j] - rho * sim[d][j];
        double fxr = objective(xr, c); fe++;
        int shrink = 0;
        if (fxr < fs[0]) {
            for (j = 0; j < d; j++) xe[j] = (1.0 + rho * chi) * xbar[j] - rho * chi * sim[d][j];
            double fxe = objective(xe, c); fe++;
            if (fxe < fxr) { memcpy(sim[d], xe, sizeof xe); fs[d] = fxe; }
            else { memcpy(sim[d], xr, sizeof xr); fs[d] = fxr; }
        } else if (fxr < fs[d - 1]) {
            memcpy(sim[d], xr, sizeof xr); fs[d] = fxr;
        } else if (fxr < fs[d]) {
            for (j = 0; j < d; j++) xc[j] = (1.0 + psi * rho) * xbar[j] - psi * rho * sim[d][j];
            double fxc = objective(xc, c); fe++;
            if (fxc <= fxr) { memcpy(sim[d], xc, sizeof xc); fs[d] = fxc; } else shrink = 1;
        } else {
            for (j = 0; j < d; j++) xc[j] = (1.0 - psi) * xbar[j] + psi * sim[d][j];
            double fxc = objective(xc, c); fe++;
            if (fxc < fs[d]) { memcpy(sim[d], xc, sizeof xc); fs[d] = fxc; } else shrink = 1;
        }
        if (shrink) {
            for (k = 1; k <= d; k++) {
                for (j = 0; j < d; j++) sim[k][j] = sim[0][j] + sig * (sim[k][j] - sim[0][j]);
                fs[k] = objective(sim[k], c); fe++;
            }
        }
        for (k = 1; k <= d; k++) {
            double fv = fs[k]; memcpy(tmp, sim[k], sizeof tmp);
            int m = k - 1;
            while (m >= 0 && fs[m] > fv) { fs[m + 1] = fs[m]; memcpy(sim[m + 1], sim[m], sizeof tmp); m--; }
            fs[m + 1] = fv; memcpy(sim[m + 1], tmp, sizeof tmp);
        }
        it++;
    }
    for (j = 0; j < d; j++) x[j] = sim[0][j];
    *nfev += fe;
    return fs[0];
}

/* Nelder-Mead and one restart from its result */
static double nm2(const fctx *c, double *x, int *nfev)
{
    int d = ctx_npar(c);
    nelder_mead(c, x, d, 8000, 1e-8, 1e-10, nfev);
    return nelder_mead(c, x, d, 8000, 1e-8, 1e-10, nfev);
}

/* a loose Nelder-Mead pass, used to screen starting values and candidate locations before the
 * full search is run from the most promising of them */
static double nm_loose(const fctx *c, double *x, int *nfev)
{
    return nelder_mead(c, x, ctx_npar(c), 1500, 1e-4, 1e-6, nfev);
}

/* indices of the (at most) nkeep smallest values among v[0..m-1] within slack of the smallest */
static int screen(const double *v, int m, double slack, int nkeep, int *keep)
{
    int nk = 0;
    double vmin = INFINITY;
    for (int i = 0; i < m; i++) if (v[i] < vmin) vmin = v[i];
    for (int r = 0; r < nkeep; r++) {
        int bi = -1;
        for (int i = 0; i < m; i++) {
            int used = 0;
            for (int q = 0; q < nk; q++) if (keep[q] == i) used = 1;
            if (used || !(v[i] <= vmin + slack)) continue;
            if (bi < 0 || v[i] < v[bi]) bi = i;
        }
        if (bi < 0) break;
        keep[nk++] = bi;
    }
    return nk;
}

/* splitmix64, for the fixed-seed random starting values of grouped fits */
static double sm64_unif(unsigned long long *s)
{
    unsigned long long z = (*s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    z ^= z >> 31;
    return ((double)(z >> 11) + 0.5) / 9007199254740992.0;
}

static double sm64_norm(unsigned long long *s)
{
    double u1 = sm64_unif(s), u2 = sm64_unif(s);
    return sqrt(-2.0 * log(u1)) * cos(TWO_PI * u2);
}

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

/* One branch of the fitting procedure. Exact data: Nelder-Mead from a grid of 27 starts about the
 * circular mean (9 and a spike on the most-tied observation at a fixed shape), under the shape
 * floor; when the fitted shape is at most one the location is profiled over the stored observations
 * within `window` of the fit (the log-likelihood is then convex in mu between observations, so every
 * maximising location is an observation), followed by a free polish. Grouped data: 30 fixed-seed
 * random starts, every fifth a spike on the modal cell. Every start is first run through a loose
 * Nelder-Mead pass and the full search (two tight passes) is run from the best few, those within a
 * log-likelihood of 2 of the best loose value; candidate locations are ranked by the likelihood at
 * the current (sigma, lam, beta), the best 12 and every tied value get a loose pass, and the full
 * search is run from those within 1 of the best. out[0..7] = mu, sigma, lam, beta, log-likelihood,
 * function evaluations, 0 (the flags are set by fit_core), number of starts. */
static void fit_branch(const fctx *c0, double window, const double *start, unsigned long long seed0,
                       int exhaustive, double *out)
{
    fctx c = *c0;
    c.fixmu_on = 0;
    int nfev = 0, nst = 0, d = ctx_npar(&c), keep[6], nk;
    (void)exhaustive;   /* kept for the interface; the screen below decides */
    double xs[40][4], vs[40], best[4] = {0.0, 0.0, 0.0, 0.0}, bestv = INFINITY, x[4];
    double bmu = 0.0, bla = 0.0, blb = 0.0, blbe = 0.0;
    if (start) {
        /* a bootstrap refit: the original fit (mu, sigma, lam, beta, log a, log b) and three random
         * perturbations of it, sd (0.05, 0.3, 0.3, 0.2) in (mu, log a, log b, log beta) */
        unsigned long long seed = seed0;
        double smu = start[0], sla = start[4], slb = start[5];
        double slbe = log(start[3]);
        for (int i = 0; i < 4; i++) {
            double e0 = 0.0, e1 = 0.0, e2 = 0.0, e3 = 0.0;
            if (i > 0) {
                e0 = 0.05 * sm64_norm(&seed); e1 = 0.3 * sm64_norm(&seed);
                e2 = 0.3 * sm64_norm(&seed); e3 = 0.2 * sm64_norm(&seed);
            }
            double lbe = slbe + e3;
            if (c.fixbeta_on) lbe = log(c.fixbeta);
            else if (exp(lbe) < c.beta0) lbe = log(1.25 * c.beta0);
            pack(xs[nst], &c, smu + e0, sla + e1, c.sym ? sla + e1 : slb + e2, lbe);
            nst++;
        }
    } else if (!c.grouped) {
        double C = 0.0, S = 0.0;
        for (int i = 0; i < c.n; i++) { C += cos(c.th[i]); S += sin(c.th[i]); }
        double m = atan2(S, C), R = sqrt(C * C + S * S) / (double)c.n;
        if (R < 1e-3) R = 1e-3;
        double s0 = sqrt(-2.0 * log(R)) / sqrt(2.0);
        if (s0 < 1e-3) s0 = 1e-3;
        const double dms[3] = {-0.3, 0.0, 0.3}, l0s[3] = {0.3, 0.5, 0.7}, b0s[3] = {0.6, 1.0, 2.0};
        double bfix = c.fixbeta_on ? c.fixbeta : 0.0;
        for (int i = 0; i < 3; i++)
            for (int l = 0; l < 3; l++) {
                if (c.sym && l > 0) continue;
                double l0 = c.sym ? 0.5 : l0s[l];
                for (int bi = 0; bi < 3; bi++) {
                    double b0 = b0s[bi];
                    if (bfix > 0.0) { if (bi > 0) continue; b0 = bfix; }
                    else if (b0 < c.beta0) b0 = 1.25 * c.beta0;
                    pack(xs[nst], &c, m + dms[i], log(l0 / s0), log((1.0 - l0) / s0), log(b0));
                    nst++;
                }
            }
        if (bfix > 0.0) {
            /* at a fixed or floor shape, also a spike on the most-tied observation, with the rates at which
             * the likelihood of a spike grows fastest, log(s / beta) / beta, s = n / (n - m) */
            double *u = (double *)malloc(sizeof(double) * (size_t)c.n);
            if (u) {
                memcpy(u, c.th, sizeof(double) * (size_t)c.n);
                qsort(u, (size_t)c.n, sizeof(double), cmp_double);
                int run = 1, bestrun = 1;
                double vt = u[0];
                for (int i = 1; i < c.n; i++) {
                    run = (u[i] == u[i - 1]) ? run + 1 : 1;
                    if (run > bestrun) { bestrun = run; vt = u[i]; }
                }
                free(u);
                double sst = (double)c.n / (double)(c.n - bestrun > 0 ? c.n - bestrun : 1);
                double lsp = log(sst / bfix) / bfix;
                if (lsp < 1.0) lsp = 1.0;
                if (lsp > 0.9 * c.rlim) lsp = 0.9 * c.rlim;
                pack(xs[nst], &c, vt, lsp, lsp, log(bfix));
                nst++;
            }
        }
        if (!c.sym) {
            /* one-sided starts: the density supported on one side of mu only (one log-rate near the
             * bound), with mu at either end of the widest gap between observations */
            double *u = (double *)malloc(sizeof(double) * (size_t)c.n);
            if (u) {
                memcpy(u, c.th, sizeof(double) * (size_t)c.n);
                qsort(u, (size_t)c.n, sizeof(double), cmp_double);
                int g = c.n - 1;
                double gmax = u[0] + TWO_PI - u[c.n - 1];
                for (int i = 0; i + 1 < c.n; i++)
                    if (u[i + 1] - u[i] > gmax) { gmax = u[i + 1] - u[i]; g = i; }
                double after = u[(g + 1) % c.n], before = u[g];
                double spread = TWO_PI - gmax;
                double lmod = -log(spread > 1e-3 ? spread : 1e-3), lbig = 0.9 * c.rlim;
                for (int bi = 0; bi < 3; bi++) {
                    double b0 = b0s[bi];
                    if (bfix > 0.0) { if (bi > 0) continue; b0 = bfix; }
                    else if (b0 < c.beta0) b0 = 1.25 * c.beta0;
                    pack(xs[nst], &c, after, lmod, lbig, log(b0));
                    nst++;
                    pack(xs[nst], &c, before, lbig, lmod, log(b0));
                    nst++;
                }
                free(u);
            }
        }
    } else {
        unsigned long long seed = 3ULL;
        int imode = 0;
        double ntot = 0.0;
        for (int j = 0; j < c.J; j++) { ntot += c.cnt[j]; if (c.cnt[j] > c.cnt[imode]) imode = j; }
        double vmode = c.orig + c.hc * (double)imode;
        for (int i = 0; i < 30; i++) {
            double lbe = log(0.05) + (log(5.0) - log(0.05)) * sm64_unif(&seed);
            double mu0 = vmode + 0.1 * sm64_norm(&seed);
            double la0 = 1.0 + sm64_norm(&seed), lb0 = 1.0 + sm64_norm(&seed);
            double be = exp(lbe);
            if (i % 5 == 0) {
                double sst = ntot / (ntot - c.cnt[imode] > 0.0 ? ntot - c.cnt[imode] : 1.0);
                mu0 = vmode; la0 = log(sst / be) / be; lb0 = la0;
            }
            if (c.fixbeta_on) lbe = log(c.fixbeta);
            else if (be < c.beta0) lbe = log(1.25 * c.beta0);
            pack(xs[nst], &c, mu0, la0, lb0, lbe);
            nst++;
        }
    }
    {
        int nfa[40];
        WSGN_TASKLOOP(c.par)
        for (int i = 0; i < nst; i++) { int nf = 0; vs[i] = nm_loose(&c, xs[i], &nf); nfa[i] = nf; }
        for (int i = 0; i < nst; i++) nfev += nfa[i];
    }
    nk = screen(vs, nst, 2.0, 3, keep);
    {
        double xq[6][4], vq[6];
        int nfq[6];
        WSGN_TASKLOOP(c.par)
        for (int q = 0; q < nk; q++) {
            int nf = 0;
            memcpy(xq[q], xs[keep[q]], sizeof xq[q]);
            vq[q] = nm2(&c, xq[q], &nf);
            nfq[q] = nf;
        }
        for (int q = 0; q < nk; q++) {
            nfev += nfq[q];
            if (vq[q] < bestv) { bestv = vq[q]; memcpy(best, xq[q], sizeof best); }
        }
    }
    unpack(best, &c, &bmu, &bla, &blb, &blbe);
    if (!c.grouped && exp(blbe) <= 1.0 && window > 0.0 && c.n > 0) {
        double *u = (double *)malloc(sizeof(double) * (size_t)c.n);
        double *cv = (double *)malloc(sizeof(double) * (size_t)c.n);
        double *cx = (double *)malloc(sizeof(double) * 4 * (size_t)c.n);
        double *cm = (double *)malloc(sizeof(double) * (size_t)c.n);
        if (u && cv && cx && cm) {
            memcpy(u, c.th, sizeof(double) * (size_t)c.n);
            qsort(u, (size_t)c.n, sizeof(double), cmp_double);
            fctx cf = c;
            cf.fixmu_on = 1;
            int nc = 0;
            int *tie = (int *)malloc(sizeof(int) * (size_t)c.n);
            double *g1 = (double *)malloc(sizeof(double) * (size_t)c.n);
            double *g2 = (double *)malloc(sizeof(double) * (size_t)c.n);
            double *c0v = (double *)malloc(sizeof(double) * (size_t)c.n);
            if (!tie || !g1 || !g2 || !c0v) nc = -1;
            /* pass 1: the distinct observations within the window, each with its likelihood at the
             * current fit and the gaps to its two neighbouring distinct values */
            for (int i = 0; i < c.n && nc >= 0; i++) {
                if (i > 0 && u[i] == u[i - 1]) continue;
                double dd = fabs(remainder(u[i] - bmu, TWO_PI));
                if (!(dd < window)) continue;
                int ip = i, inx = i;
                do { ip = (ip + c.n - 1) % c.n; } while (u[ip] == u[i] && ip != i);
                do { inx = (inx + 1) % c.n; } while (u[inx] == u[i] && inx != i);
                double dp = fabs(remainder(u[i] - u[ip], TWO_PI)), dx = fabs(remainder(u[inx] - u[i], TWO_PI));
                g1[nc] = (dp < dx) ? dp : dx;
                g2[nc] = (dp < dx) ? dx : dp;
                cf.fixmu = u[i];
                pack(cx + 4 * nc, &cf, u[i], bla, blb, blbe);
                c0v[nc] = objective(cx + 4 * nc, &cf);
                nfev++;
                cm[nc] = u[i];
                tie[nc] = (i + 1 < c.n && u[i + 1] == u[i]);
                nc++;
            }
            /* pass 2: the candidates to screen: all of them when there are at most 200, otherwise
             * the 15 best at the current fit, the 15 with the smallest gap to a neighbour, the 15
             * with the smallest gap to the further neighbour, and every tied value */
            int *sel = (nc > 0) ? (int *)malloc(sizeof(int) * (size_t)nc) : NULL;
            int ns = 0;
            if (sel) {
                if (nc <= 200) { for (int q = 0; q < nc; q++) sel[ns++] = q; }
                else {
                    char *mk = (char *)calloc((size_t)nc, 1);
                    if (mk) {
                        const double *crit[3] = {c0v, g1, g2};
                        for (int cr = 0; cr < 3; cr++) {
                            for (int r = 0; r < 15; r++) {
                                int bi = -1;
                                for (int q = 0; q < nc; q++)
                                    if (!mk[q] && (bi < 0 || crit[cr][q] < crit[cr][bi])) bi = q;
                                if (bi < 0) break;
                                mk[bi] = 1;
                            }
                        }
                        for (int q = 0; q < nc; q++) if (mk[q] || tie[q]) sel[ns++] = q;
                        free(mk);
                    } else { for (int q = 0; q < nc; q++) sel[ns++] = q; }
                }
            }
            /* screen each selected candidate after a short optimisation of the other parameters,
             * from the current fit and from a spike as wide as the gap to the nearest other value
             * (which finds spikes on near-ties), keeping the better */
            int *nfs = (ns > 0) ? (int *)calloc((size_t)ns, sizeof(int)) : NULL;
            if (ns > 0 && !nfs) ns = 0;
            WSGN_TASKLOOP(c.par)
            for (int q0 = 0; q0 < ns; q0++) {
                int q = sel[q0], dq, nf = 0;
                fctx cl = cf;
                double xa[4], lsp = -log(g1[q] > 1e-300 ? g1[q] : 1.0);
                if (lsp > 0.9 * c.rlim) lsp = 0.9 * c.rlim;
                if (lsp < -0.9 * c.rlim) lsp = -0.9 * c.rlim;
                cl.fixmu = cm[q];
                dq = ctx_npar(&cl);
                cv[q] = nelder_mead(&cl, cx + 4 * q, dq, 20 * dq, 1e-3, 1e-4, &nf);
                pack(xa, &cl, cm[q], lsp, lsp, blbe);
                double va = nelder_mead(&cl, xa, dq, 20 * dq, 1e-3, 1e-4, &nf);
                if (va < cv[q]) { cv[q] = va; memcpy(cx + 4 * q, xa, sizeof xa); }
                nfs[q0] = nf;
            }
            for (int q0 = 0; q0 < ns; q0++) nfev += nfs[q0];
            free(nfs);
            if (nc > 0 && ns > 0) {
                int kk[6], n12 = 0;
                int *k12 = (int *)malloc(sizeof(int) * (size_t)ns);
                double *v12 = (double *)malloc(sizeof(double) * (size_t)ns);
                double *vs2 = (double *)malloc(sizeof(double) * (size_t)ns);
                if (!k12 || !v12 || !vs2) { free(k12); free(v12); free(vs2); k12 = NULL; v12 = NULL; vs2 = NULL; ns = 0; }
                if (ns > 0) {
                    /* the 12 best screened candidates and every tied value get a loose pass */
                    for (int q0 = 0; q0 < ns; q0++) vs2[q0] = cv[sel[q0]];
                    n12 = screen(vs2, ns, INFINITY, 12, k12);
                    for (int q0 = 0; q0 < ns; q0++) {
                        if (!tie[sel[q0]]) continue;
                        int seen = 0;
                        for (int r2 = 0; r2 < n12; r2++) if (k12[r2] == q0) seen = 1;
                        if (!seen) k12[n12++] = q0;
                    }
                    for (int q = 0; q < n12; q++) k12[q] = sel[k12[q]];
                    int *nf12 = (int *)calloc((size_t)(n12 > 0 ? n12 : 1), sizeof(int));
                    if (!nf12) n12 = 0;
                    WSGN_TASKLOOP(c.par)
                    for (int q = 0; q < n12; q++) {
                        int nf = 0;
                        fctx cl = cf;
                        cl.fixmu = cm[k12[q]];
                        v12[q] = nm_loose(&cl, cx + 4 * k12[q], &nf);
                        nf12[q] = nf;
                    }
                    for (int q = 0; q < n12; q++) nfev += nf12[q];
                    free(nf12);
                }
                free(vs2);
                int nkk = (nc > 0) ? screen(v12, n12, 1.0, 2, kk) : 0;
                for (int q = 0; q < nkk; q++) kk[q] = k12[kk[q]];
                free(k12); free(v12);
                {
                    double xk[2][4], vk[2];
                    int nfk[2];
                    WSGN_TASKLOOP(c.par)
                    for (int q = 0; q < nkk; q++) {
                        int nf = 0;
                        fctx cl = cf;
                        cl.fixmu = cm[kk[q]];
                        memcpy(xk[q], cx + 4 * kk[q], sizeof xk[q]);
                        vk[q] = nm2(&cl, xk[q], &nf);
                        nfk[q] = nf;
                    }
                    for (int q = 0; q < nkk; q++) {
                        nfev += nfk[q];
                        if (vk[q] < bestv) {
                            double mm;
                            bestv = vk[q];
                            unpack(xk[q], &cf, &mm, &bla, &blb, &blbe);
                            bmu = cm[kk[q]];
                        }
                    }
                }
            }
            free(tie); free(g1); free(g2); free(c0v); free(sel);
            pack(x, &c, bmu, bla, blb, blbe);
            double v = nm2(&c, x, &nfev);
            if (v < bestv - 1e-9) { bestv = v; unpack(x, &c, &bmu, &bla, &blb, &blbe); }
        }
        free(u); free(cv); free(cx); free(cm);
    }
    (void)d;
    double a = exp(bla), b = exp(blb), be = exp(blbe);
    out[0] = mod2pi(bmu);
    out[1] = 1.0 / (a + b);
    out[2] = 1.0 / (1.0 + exp(blb - bla));
    out[3] = be;
    out[4] = -bestv;
    out[5] = (double)nfev;
    out[6] = 0.0;
    out[7] = (double)nst;
    out[10] = bla;
    out[11] = blb;
}

/* The fit: the procedure above, and for exact data with the shape free also a search of the spike
 * branch, whose likelihood grows as the shape falls and whose profile in the shape can be ragged
 * just above the floor: the same procedure at the fixed shapes beta0 (1, 1.25, 1.6), capped
 * below beta1/2, each with a spike start on the most-tied observation and the location searched over
 * every observation, followed by a free polish (shape above the floor) from each of them. The best of everything is returned, so the result is
 * the maximum over [beta0, beta1] found by the regular and the spike searches, whether the spike
 * maximum lies on the floor or a little above it. out[0..9] = mu, sigma, lam, beta, log-likelihood,
 * evaluations, at-floor flag, starts, at-upper-bound flag, at-rate-bound flag. */
static void fit_core(const fctx *c0, double window, const double *start, unsigned long long seed0,
                     int exhaustive, double *out)
{
    if (!c0->grouped && !c0->fixbeta_on) {
        const double fac[3] = {1.0, 1.25, 1.6};
        double og[3][NOUT];
        int ng = 0;
        for (int g = 0; g < 3; g++) {
            double bg = c0->beta0 * fac[g];
            if (g > 0 && bg >= 0.5 * c0->beta1) break;
            ng++;
        }
        /* the regular branch and the spike branches are independent of one another, and the free
         * polish of each spike fit depends on that fit alone, so it follows it in the same task */
        double xp[3][4], vp[3];
        int nfp[3] = {0, 0, 0};
        WSGN_TASKLOOP(c0->par)
        for (int g = 0; g <= ng; g++) {
            if (g == 0) fit_branch(c0, window, start, seed0, exhaustive, out);
            else {
                fctx cf = *c0;
                cf.fixbeta_on = 1;
                cf.fixbeta = c0->beta0 * fac[g - 1];
                double *o2 = og[g - 1];
                fit_branch(&cf, M_PI, NULL, 0ULL, 0, o2);
                if (o2[4] > -INFINITY) {
                    fctx cp = *c0;
                    double mu = o2[0], la = log(o2[2] / o2[1]), lb = log((1.0 - o2[2]) / o2[1]);
                    double lbe = log(o2[3]);
                    int nf = 0;
                    cp.fixmu_on = 0;
                    pack(xp[g - 1], &cp, mu, la, lb, lbe);
                    vp[g - 1] = nm2(&cp, xp[g - 1], &nf);
                    nfp[g - 1] = nf;
                }
            }
        }
        for (int g = 0; g < ng; g++) {
            out[5] += og[g][5];
            out[7] += og[g][7];
        }
        /* free polish from every grid fit, best first */
        int order[3] = {0, 1, 2};
        for (int i = 0; i < ng; i++)
            for (int j = i + 1; j < ng; j++)
                if (og[order[j]][4] > og[order[i]][4]) { int t = order[i]; order[i] = order[j]; order[j] = t; }
        for (int q = 0; q < ng; q++) {
            int g = order[q];
            double *o2 = og[g];
            if (!(o2[4] > -INFINITY)) continue;
            fctx cp = *c0;
            double *x = xp[g], v = vp[g], mu, la, lb, lbe;
            cp.fixmu_on = 0;
            out[5] += (double)nfp[g];
            if (o2[4] > out[4]) {
                for (int j = 0; j < 5; j++) out[j] = o2[j];
                out[10] = o2[10]; out[11] = o2[11];
            }
            if (-v > out[4] + 1e-9) {
                unpack(x, &cp, &mu, &la, &lb, &lbe);
                out[0] = mod2pi(mu);
                out[1] = 1.0 / (exp(la) + exp(lb));
                out[2] = 1.0 / (1.0 + exp(lb - la));
                out[3] = exp(lbe);
                out[4] = -v;
                out[10] = la; out[11] = lb;
            }
        }
    } else {
        fit_branch(c0, window, start, seed0, exhaustive, out);
    }
    /* A fit that ends on the upper bound of the shape may have stopped in a two-sided local maximum
     * there while the maximum at the bound is one-sided; the search at the bound itself (the
     * profile at beta1, with its one-sided starts and location screen) is kept if higher. */
    if (!c0->fixbeta_on && isfinite(out[4]) && out[3] >= c0->beta1 * (1.0 - 1e-5)) {
        fctx cb = *c0;
        double ob[NOUT];
        cb.fixbeta_on = 1;
        cb.fixbeta = c0->beta1;
        fit_branch(&cb, window, NULL, 0ULL, 0, ob);
        out[5] += ob[5];
        out[7] += ob[7];
        if (ob[4] > out[4] + 1e-9) {
            for (int j = 0; j < 5; j++) out[j] = ob[j];
            out[10] = ob[10]; out[11] = ob[11];
        }
    }
    /* A one-sided fit (one log-rate effectively infinite) sits on a cusp in the location, where any
     * step of mu loses by orders of magnitude, so the free polishes above can collapse their simplex
     * along mu before the other parameters have moved. One more polish with mu fixed at the fitted
     * location and the larger log-rate pulled back to 45 above the smaller; kept only if higher. */
    if (fabs(out[10] - out[11]) > 27.631021115928547 && isfinite(out[4])) {
        fctx cp = *c0;
        double x[4], mu = out[0], la = out[10], lb = out[11], lbe = log(out[3]);
        double lsm = (la < lb) ? la : lb, lnew = lsm + 45.0;
        if (la > lb) { if (lnew < la) la = lnew; else lnew = la; }
        else { if (lnew < lb) lb = lnew; else lnew = lb; }
        if (fabs(lnew) <= c0->rlim) {
            int nf = 0;
            cp.fixmu_on = 1;
            cp.fixmu = mu;
            pack(x, &cp, mu, la, lb, lbe);
            double v = nm2(&cp, x, &nf);
            out[5] += (double)nf;
            if (-v > out[4] + 1e-9) {
                unpack(x, &cp, &mu, &la, &lb, &lbe);
                out[1] = 1.0 / (exp(la) + exp(lb));
                out[2] = 1.0 / (1.0 + exp(lb - la));
                out[3] = exp(lbe);
                out[4] = -v;
                out[10] = la; out[11] = lb;
            }
        }
    }
    double la = out[10], lb = out[11];
    /* at a shape up to one every maximising location is an observation; a polish can drift a few
     * ulps off it through rounding, so snap back when that costs nothing */
    if (!c0->grouped && out[3] <= 1.0) {
        int inear = -1;
        double dnear = 1e-9;
        for (int i = 0; i < c0->n; i++) {
            double dd = fabs(remainder(c0->th[i] - out[0], TWO_PI));
            if (dd < dnear) { dnear = dd; inear = i; }
        }
        if (inear >= 0 && c0->th[inear] != out[0]) {
            double v = loglik_ctx(c0, c0->th[inear], la, lb, out[3]);
            if (v >= out[4] - 1e-9) { out[0] = c0->th[inear]; out[4] = v; }
        }
    }
    out[6] = (!c0->fixbeta_on && out[3] <= c0->beta0 * (1.0 + 1e-6)) ? 1.0 : 0.0;
    out[8] = (!c0->fixbeta_on && out[3] >= c0->beta1 * (1.0 - 1e-5)) ? 1.0 : 0.0;
    out[9] = (fmax(fabs(la), fabs(lb)) >= 0.999 * c0->rlim) ? 1.0 : 0.0;
}

/* bin angles into the J cells of width hc centred at orig + j hc; non-finite angles are skipped
 * (the front ends refuse them), and the index is reduced modulo J before the cast, so it is always
 * in range */
static void bin_cells(const double *th, int n, double hc, double orig, int J, double *cnt)
{
    for (int j = 0; j < J; j++) cnt[j] = 0.0;
    for (int i = 0; i < n; i++) {
        double z = floor((th[i] - orig) / hc + 0.5);
        if (!isfinite(z)) continue;
        double r = fmod(z, (double)J);
        if (r < 0.0) r += (double)J;
        long jj = (long)r;
        if (jj < 0 || jj >= J) continue;
        cnt[jj] += 1.0;
    }
}

#define WSGN_MAXCELLS 1e6

static void setup_ctx(fctx *c, const double *th, int n, int grouped, double res, double orig,
                      double beta0, int sym, double fixbeta, int K)
{
    memset(c, 0, sizeof *c);
    c->th = th; c->n = n; c->grouped = grouped; c->sym = sym; c->K = K;
    c->hc = res; c->orig = orig;
    /* more than WSGN_MAXCELLS cells (or a non-finite count) is refused: J stays 0 and the entry
     * points return NaN, so the cast below is never out of range */
    {
        double Jd = grouped ? floor(TWO_PI / res + 0.5) : 0.0;
        c->J = (Jd >= 1.0 && Jd <= WSGN_MAXCELLS) ? (int)Jd : 0;
    }
    c->beta0 = beta0;
    c->beta1 = grouped ? 50.0 : 20.0;
    /* bound on |log rate|: 60 for grouped data; for exact data large enough for the spike at the
     * smallest shape in play (its rates are near log(s/beta)/beta) and at least 30 */
    {
        /* the floor sets the bound; a fixed shape below the floor (a profile there) widens it */
        double bb = (fixbeta > 0.0 && fixbeta < beta0) ? fixbeta : beta0;
        double r = 3.0 * (1.0 - log(bb)) / bb;
        c->rlim = grouped ? 60.0 : (r > 30.0 ? r : 30.0);
    }
    c->fixbeta_on = (fixbeta > 0.0);
    c->fixbeta = fixbeta;
}

static void fit_one(const double *th, int n, int grouped, double res, double orig, double beta0,
                    int sym, double fixbeta, double window, int K, const double *start,
                    unsigned long long seed, int exhaustive, int par, double *out)
{
    fctx c;
    double *tr = (double *)malloc(sizeof(double) * (size_t)(n > 0 ? n : 1));
    if (!tr) { for (int i = 0; i < NOUT; i++) out[i] = NAN; return; }
    for (int i = 0; i < n; i++) tr[i] = mod2pi(th[i]);
    setup_ctx(&c, tr, n, grouped, res, orig, beta0, sym, fixbeta, K);
    c.par = par;
    double *cnt = NULL;
    if (grouped && c.J < 1) { free(tr); for (int i = 0; i < NOUT; i++) out[i] = NAN; return; }
    if (grouped) {
        cnt = (double *)malloc(sizeof(double) * (size_t)c.J);
        if (!cnt) { free(tr); for (int i = 0; i < NOUT; i++) out[i] = NAN; return; }
        bin_cells(tr, n, res, orig, c.J, cnt);
        c.cnt = cnt;
    }
    fit_core(&c, window, start, seed, exhaustive, out);
    if (cnt) free(cnt);
    free(tr);
}

/* The independent trials inside a fit (the starts, the candidate locations, the spike branches)
 * run as OpenMP tasks on `*nthreads` threads when the library is built with OpenMP; the best is
 * chosen afterwards in index order, so the result is the same for every number of threads. */
void wsgn_fit(const double *theta, const int *n, const int *grouped, const double *res,
              const double *origin, const double *beta0, const int *sym, const double *fixbeta,
              const double *window, const int *K, const int *nthreads, double *out)
{
#ifdef WSGN_TASKS
    int nth = (*nthreads > 1) ? *nthreads : 1;
    if (nth > 1 && !omp_in_parallel()) {
#pragma omp parallel num_threads(nth)
#pragma omp single
        fit_one(theta, *n, *grouped, *res, *origin, *beta0, *sym, *fixbeta, *window, *K, NULL,
                0ULL, 0, 1, out);
        return;
    }
#else
    (void)nthreads;
#endif
    fit_one(theta, *n, *grouped, *res, *origin, *beta0, *sym, *fixbeta, *window, *K, NULL, 0ULL, 0,
            0, out);
}

/* beta0 is the shape floor of the fit to compare with: it sets the rate bound (a shape below it
 * widens the bound, as in setup_ctx), whatever the shape being profiled */
void wsgn_profile(const double *theta, const int *n, const int *grouped, const double *res,
                  const double *origin, const double *beta0, const double *betas, const int *nb,
                  const int *sym, const double *window, const int *K, const int *nthreads,
                  double *out)
{
    int nbeta = *nb;
#ifdef WSGN_TASKS
    int nth = (*nthreads > 1) ? *nthreads : 1;
    if (nth > 1 && !omp_in_parallel()) {
#pragma omp parallel num_threads(nth)
#pragma omp single
        {
            /* the shapes are independent of one another, and each fit shares the threads */
            WSGN_TASKLOOP(1)
            for (int i = 0; i < nbeta; i++)
                fit_one(theta, *n, *grouped, *res, *origin, *beta0, *sym, betas[i], *window, *K,
                        NULL, 0ULL, 0, 1, out + NOUT * i);
        }
        return;
    }
#else
    (void)nthreads;
#endif
    for (int i = 0; i < nbeta; i++)
        fit_one(theta, *n, *grouped, *res, *origin, *beta0, *sym, betas[i], *window, *K, NULL,
                0ULL, 0, 0, out + NOUT * i);
}

void wsgn_grouped_ll(const double *theta, const int *n, const double *res, const double *origin,
                     const double *mu, const double *sigma, const double *beta, const double *lam,
                     const int *K, double *out)
{
    fctx c;
    if (!(*res > 0.0 && *res <= M_PI) || *n < 0) { *out = NAN; return; }
    setup_ctx(&c, theta, *n, 1, *res, *origin, 0.0, 0, 0.0, *K);
    if (c.J < 1) { *out = NAN; return; }
    double *cnt = (double *)malloc(sizeof(double) * (size_t)c.J);
    if (!cnt) { *out = NAN; return; }
    bin_cells(theta, *n, *res, *origin, c.J, cnt);
    c.cnt = cnt;
    *out = loglik_ctx(&c, *mu, log(*lam / *sigma), log((1.0 - *lam) / *sigma), *beta);
    free(cnt);
}

/* Parametric-bootstrap refits: column b of samples (n x B, column-major) is refitted by the same
 * procedure, started from the fit par0 and three fixed-seed random
 * perturbations of it, with the same location search and floor branch; est (10 x B) receives the
 * fits, se[0..3] the bootstrap standard deviations of (mu, sigma, lam, beta), mu as circular deviations from par0[0], and se[4] the fraction of refits
 * with shape below one. Replicates are independent, so they run in parallel under OpenMP when the
 * library is built with it, with results that do not depend on the number of threads. */
void wsgn_fit_batch(const double *samples, const int *n, const int *B, const int *grouped,
                    const double *res, const double *origin, const double *beta0, const int *sym,
                    const double *window, const int *K, const int *nthreads, const double *par0,
                    double *est, double *se)
{
    int nb = *B;
#ifdef _OPENMP
    int nth = (*nthreads > 0) ? *nthreads : 1;
#pragma omp parallel for schedule(dynamic) num_threads(nth)
#endif
    for (int b = 0; b < nb; b++)
        fit_one(samples + (size_t)b * (size_t)(*n), *n, *grouped, *res, *origin, *beta0, *sym, 0.0,
                *window, *K, par0, 5000ULL + (unsigned long long)b, 0, 0,
                est + NOUT * (size_t)b);
    (void)nthreads;
    for (int j = 0; j < 4; j++) {
        double m = 0.0, v = 0.0;
        for (int b = 0; b < nb; b++) {
            double e = est[NOUT * b + j];
            if (j == 0) e = remainder(e - par0[0], TWO_PI);
            m += e;
        }
        m /= (double)nb;
        for (int b = 0; b < nb; b++) {
            double e = est[NOUT * b + j];
            if (j == 0) e = remainder(e - par0[0], TWO_PI);
            v += (e - m) * (e - m);
        }
        se[j] = (nb > 1) ? sqrt(v / (double)(nb - 1)) : NAN;
    }
    double f = 0.0;
    for (int b = 0; b < nb; b++) if (est[NOUT * b + 3] < 1.0) f += 1.0;
    se[4] = f / (double)nb;
}

/* Expected Fisher information per observation of the wrapped distribution in (mu, sigma, lam,
 * beta), I = int s s' g, with the scores s by central differences on steps relative to each
 * parameter's own scale (1e-5 of sigma for mu and sigma, of min(lam, 1 - lam) for lam, of beta for
 * beta), so that a concentrated or nearly one-sided fit keeps its accuracy. The integral is the
 * midpoint rule in v on theta = mu + c sinh(v), |theta - mu| < pi, with c = min(1, 0.5 / max(a, b)),
 * which puts most of the ngrid points within a few scales of mu however small sigma is, and never
 * a point at mu itself (ngrid even). vcov = I^-1 / n and se its root diagonal. For the symmetric
 * member (lam = 1/2 fixed) the lam row and column are NaN. The information is finite and
 * nonsingular for beta > 1/2; the normal approximation holds for beta > 1. */
void wsgn_info(const double *par, const int *sym, const int *K, const int *ngrid, const int *n,
               double *info, double *vcov, double *se)
{
    int act[4], na = 0;
    for (int j = 0; j < 4; j++) if (!(*sym && j == 2)) act[na++] = j;
    double I[4][4] = {{0}}, p0[4], hstep[4];
    for (int j = 0; j < 4; j++) p0[j] = par[j];
    double lmin = (p0[2] < 1.0 - p0[2]) ? p0[2] : 1.0 - p0[2];
    hstep[0] = 1e-5 * (p0[1] < 1.0 ? p0[1] : 1.0);
    hstep[1] = 1e-5 * p0[1];
    hstep[2] = 1e-5 * lmin;
    hstep[3] = 1e-5 * p0[3];
    for (int j = 0; j < 16; j++) { info[j] = NAN; vcov[j] = NAN; }
    for (int j = 0; j < 4; j++) se[j] = NAN;
    if (!(p0[1] > 0.0 && lmin > 0.0 && p0[3] > 0.0 && isfinite(p0[0]) && isfinite(p0[1]) &&
          isfinite(p0[3]))) return;
    double a = p0[2] / p0[1], b = (1.0 - p0[2]) / p0[1];
    double cc = 0.5 / (a > b ? a : b);
    if (cc > 1.0) cc = 1.0;
    int N = *ngrid + (*ngrid % 2);
    double V = asinh(M_PI / cc), dv = 2.0 * V / (double)N;
    for (int i = 0; i < N; i++) {
        double v = -V + ((double)i + 0.5) * dv;
        double th = p0[0] + cc * sinh(v);
        double s[4];
        double lg0 = logdens_rates(th, p0[0], log(p0[2] / p0[1]), log((1.0 - p0[2]) / p0[1]), p0[3], *K);
        for (int q = 0; q < na; q++) {
            int j = act[q];
            double pp[4], pm[4];
            memcpy(pp, p0, sizeof pp); memcpy(pm, p0, sizeof pm);
            pp[j] += hstep[j]; pm[j] -= hstep[j];
            double lp = logdens_rates(th, pp[0], log(pp[2] / pp[1]), log((1.0 - pp[2]) / pp[1]), pp[3], *K);
            double lm = logdens_rates(th, pm[0], log(pm[2] / pm[1]), log((1.0 - pm[2]) / pm[1]), pm[3], *K);
            s[q] = (lp - lm) / (2.0 * hstep[j]);
        }
        double g = exp(lg0) * cc * cosh(v) * dv;
        for (int q = 0; q < na; q++)
            for (int r = 0; r < na; r++) I[q][r] += s[q] * s[r] * g;
    }
    for (int q = 0; q < na; q++)
        for (int r = 0; r < na; r++) if (!isfinite(I[q][r])) return;
    /* invert the na x na matrix by Gauss-Jordan elimination with partial pivoting */
    double M[4][8];
    for (int a = 0; a < na; a++)
        for (int b = 0; b < 2 * na; b++) M[a][b] = (b < na) ? I[a][b] : (b - na == a ? 1.0 : 0.0);
    int ok = 1;
    for (int col = 0; col < na && ok; col++) {
        int piv = col;
        for (int r = col + 1; r < na; r++) if (fabs(M[r][col]) > fabs(M[piv][col])) piv = r;
        if (!(fabs(M[piv][col]) > 0.0)) { ok = 0; break; }
        if (piv != col) for (int b = 0; b < 2 * na; b++) { double t = M[col][b]; M[col][b] = M[piv][b]; M[piv][b] = t; }
        double dv = M[col][col];
        for (int b = 0; b < 2 * na; b++) M[col][b] /= dv;
        for (int r = 0; r < na; r++) {
            if (r == col) continue;
            double f = M[r][col];
            for (int b = 0; b < 2 * na; b++) M[r][b] -= f * M[col][b];
        }
    }
    for (int a = 0; a < na; a++)
        for (int b = 0; b < na; b++) {
            info[act[a] + 4 * act[b]] = I[a][b];
            if (ok) vcov[act[a] + 4 * act[b]] = M[a][na + b] / (double)(*n);
        }
    if (ok) for (int a = 0; a < na; a++) se[act[a]] = sqrt(vcov[5 * act[a]]);
}
