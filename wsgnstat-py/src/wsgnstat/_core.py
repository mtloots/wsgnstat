"""wsgnstat: Python (ctypes) binding to the shared C back-end for the two-piece exponential power
distribution on the line and its wrapped version on the circle, treated as a uniform scale mixture
with generalised gamma mixing (the package keeps the authors' earlier name, wrapped skew generalised
normal). Pure standard library. An installed copy carries a library built by pip; a source checkout
compiles the C source on first import. The same source backs the R package 'wsgnstat', and every
loop, optimiser and bootstrap refit runs in C."""
import ctypes as _ct
import glob as _glob
import math as _math
import os as _os
import random as _random
import sys as _sys

_here = _os.path.dirname(_os.path.abspath(__file__))


# LOADING has two paths, as in arcstat-py and exactcond-py. An INSTALLED package carries a library
# that pip built at install time, so no compiler is needed. A SOURCE CHECKOUT has the .c file but no
# library, so it is compiled on first import into a writable cache directory.
def _find_built():
    for pat in ("_libwsgnstat*.so", "_libwsgnstat*.dylib", "_libwsgnstat*.pyd", "_libwsgnstat*.dll"):
        hits = sorted(_glob.glob(_os.path.join(_here, pat)))
        if hits:
            return hits[0]
    return None


def _cache_dir():
    if _os.access(_here, _os.W_OK):
        return _here
    base = _os.environ.get("XDG_CACHE_HOME") or _os.path.join(_os.path.expanduser("~"), ".cache")
    d = _os.path.join(base, "wsgnstat")
    _os.makedirs(d, exist_ok=True)
    return d


def _ext():
    return "dll" if _sys.platform == "win32" else ("dylib" if _sys.platform == "darwin" else "so")


def _compile_from_source():
    import subprocess as _sp
    srcs = sorted(f for f in _glob.glob(_os.path.join(_here, "*.c"))
                  if _os.path.basename(f) != "init.c")
    if not srcs:
        return None
    out = _os.path.join(_cache_dir(), "libwsgnstat." + _ext())
    cc = _os.environ.get("CC") or "cc"
    base = [cc, "-O2", "-fPIC", "-shared", "-I", _here] + srcs + ["-lm", "-o", out]
    try:
        try:
            if _sys.platform.startswith("linux") and _os.environ.get("WSGN_OPENMP", "auto") != "0":
                _sp.check_call(base + ["-fopenmp"], stderr=_sp.DEVNULL)
            else:
                raise OSError("serial")
        except (OSError, _sp.CalledProcessError):
            _sp.check_call(base)
    except (OSError, _sp.CalledProcessError) as exc:
        raise ImportError(
            "wsgnstat could not build its C back-end from the bundled sources (%s). In an "
            "installed copy this should never happen, because pip builds the library at install "
            "time. In a source checkout a C compiler is required; set CC if it is not 'cc'." % exc)
    return out


_libpath = _find_built()
if _libpath is None:
    _cached = _os.path.join(_cache_dir(), "libwsgnstat." + _ext())
    _srcs = [f for f in _glob.glob(_os.path.join(_here, "*.[ch]"))]
    _stale = (not _os.path.exists(_cached) or
              any(_os.path.getmtime(f) > _os.path.getmtime(_cached) for f in _srcs))
    _libpath = _compile_from_source() if _stale else _cached
if _libpath is None:
    raise ImportError("wsgnstat: no compiled back-end and no C sources to build one from.")

_lib = _ct.CDLL(_libpath)


def _has_openmp():
    # the OpenMP runtime is a dependency of the library exactly when it was built with OpenMP
    try:
        _lib.omp_get_max_threads
    except AttributeError:
        return False
    return True


HAS_OPENMP = _has_openmp()
_dp, _ip = _ct.POINTER(_ct.c_double), _ct.POINTER(_ct.c_int)
_lib.wsgn_density.argtypes = [_dp, _ip, _dp, _dp, _dp, _dp, _dp, _dp]
_lib.wsgn_pdist.argtypes = [_dp, _ip, _dp, _dp, _dp, _dp, _dp, _dp]
_lib.wsgn_cf.argtypes = [_dp, _ip, _dp, _dp, _dp, _dp, _dp, _dp, _ip, _dp, _dp]
_lib.wsgn_trigmom.argtypes = [_ip, _dp, _dp, _dp, _dp, _dp, _dp, _dp]
_lib.wsgn_decay.argtypes = [_dp, _dp, _dp, _dp]
_lib.wsgn_wdensity.argtypes = [_dp, _ip, _dp, _dp, _dp, _dp, _ip, _ip, _dp]
_lib.wsgn_nll.argtypes = [_dp, _ip, _dp, _dp, _dp, _dp, _ip, _dp]
_lib.wsgn_wcdf.argtypes = [_dp, _ip, _dp, _dp, _dp, _dp, _ip, _dp]
_lib.wsgn_wquant.argtypes = [_dp, _ip, _dp, _dp, _dp, _dp, _ip, _dp]
_lib.wsgn_fit.argtypes = [_dp, _ip, _ip, _dp, _dp, _dp, _ip, _dp, _dp, _ip, _ip, _dp]
_lib.wsgn_profile.argtypes = [_dp, _ip, _ip, _dp, _dp, _dp, _dp, _ip, _ip, _dp, _ip, _ip,
                              _dp]
_lib.wsgn_grouped_ll.argtypes = [_dp, _ip, _dp, _dp, _dp, _dp, _dp, _dp, _ip, _dp]
_lib.wsgn_fit_batch.argtypes = [_dp, _ip, _ip, _ip, _dp, _dp, _dp, _ip, _dp, _ip, _ip, _dp,
                                _dp, _dp]
_lib.wsgn_info.argtypes = [_dp, _ip, _ip, _ip, _ip, _dp, _dp, _dp]

_TWO_PI = 2.0 * _math.pi
_NAMES = ("mu", "sigma", "lam", "beta")


def _d(seq):
    seq = [float(v) for v in seq]
    a = (_ct.c_double * max(len(seq), 1))()
    for i, v in enumerate(seq):
        a[i] = v
    return a


def _i(v):
    a = (_ct.c_int * 1)()
    a[0] = int(v)
    return a


def _d1(v):
    a = (_ct.c_double * 1)()
    a[0] = float(v)
    return a


def _seq(x):
    try:
        return list(x)
    except TypeError:
        return [x]


def _sigma(sigma, lam):
    """The scale of the reciprocal-rate subfamily (a b = 1) when sigma is not given."""
    return _math.sqrt(lam * (1.0 - lam)) if sigma is None else float(sigma)


def _check_K(K):
    try:
        k = float(K)
    except (TypeError, ValueError):
        k = float("nan")
    if not (_math.isfinite(k) and k == int(k) and 1 <= k <= 1e6):
        raise ValueError("'K' must be a whole number between 1 and 1e6")


def _check_tol(tol):
    try:
        t = float(tol)
    except (TypeError, ValueError):
        t = float("nan")
    if not (_math.isfinite(t) and 0 < t < 1):
        raise ValueError("'tol' must be one finite number in (0, 1)")
    if t < 1e-14:
        import warnings as _w
        _w.warn("'tol' below 1e-14 is beyond the rounding of the quadrature; 1e-14 is used")
        t = 1e-14
    return t


def _check(sigma, beta, lam, K=1):
    if not (0 < lam < 1):
        raise ValueError("'lam' must lie in (0, 1)")
    if not (sigma > 0 and _math.isfinite(sigma)):
        raise ValueError("'sigma' must be one positive number")
    if not (beta > 0 and _math.isfinite(beta)):
        raise ValueError("'beta' must be one positive number")
    _check_K(K)


def _check_line(sigma, beta, lam1, lam2):
    for v in (sigma, beta, lam1, lam2):
        if not (v > 0 and _math.isfinite(v)):
            raise ValueError("'sigma', 'beta', 'lam1' and 'lam2' must each be one positive finite "
                             "number")


# ------------------------------------------------------------------ the line distribution
def dsgn(x, mu=0.0, sigma=1.0, beta=2.0, lam1=0.5, lam2=0.5):
    """Density of the two-piece exponential power distribution on the line at x: right rate
    lam1/sigma, left rate lam2/sigma. With lam1 = lam and lam2 = 1 - lam it is the distribution of
    mu + Z V, Z^beta ~ Gamma(1 + 1/beta), V uniform on [-sigma/(1 - lam), sigma/lam]."""
    _check_line(sigma, beta, lam1, lam2)
    arr = _d(_seq(x))
    n = len(_seq(x))
    out = (_ct.c_double * max(n, 1))()
    _lib.wsgn_density(arr, _i(n), _d1(mu), _d1(sigma), _d1(beta), _d1(lam1), _d1(lam2), out)
    return list(out)[:n]


def psgn(q, mu=0.0, sigma=1.0, beta=2.0, lam1=0.5, lam2=0.5):
    """Distribution function of the line distribution at q, in closed form (regularised incomplete
    gamma function)."""
    _check_line(sigma, beta, lam1, lam2)
    n = len(_seq(q))
    arr = _d(_seq(q))
    out = (_ct.c_double * max(n, 1))()
    _lib.wsgn_pdist(arr, _i(n), _d1(mu), _d1(sigma), _d1(beta), _d1(lam1), _d1(lam2), out)
    return list(out)[:n]


def cfsgn(t, mu=0.0, sigma=1.0, beta=2.0, lam1=0.5, lam2=0.5, tol=1e-12, method="quadrature"):
    """Characteristic function of the line distribution at t, as a list of complex numbers.

    method "quadrature" (default): the integral G(w) = int_0^inf exp(i w u - u^beta) du on a ray in
    the upper half-plane, by the trapezoid rule in log-radius, halving the step until successive sums
    agree to tol, one finite number in (0, 1); values below 1e-14, beyond the rounding of the
    quadrature, are raised to 1e-14 with a warning. Accuracy: tol is relative to the size of each
    half, G(w) / Gamma(1 + 1/beta),
    which is of order min(1, 1/|w|); the halves cancel to leading order at large |t|, so the error
    of the characteristic function is absolute, of order tol min(1, max(lam1, lam2)/(sigma |t|)) and
    (in practice of order tol: about 1e-15 to 3e-15 for shapes up to about 3, a few times 1e-15 for shapes 5
    to 20, and up to about 1e-14 for shapes up to 50, at the default tol), and values much smaller than that carry no relative accuracy
    (for the normal member, beyond t of about 8 in standard-deviation units). Against 120-digit
    series values at t = 1, 2, 3, lam1 = 0.3, lam2 = 0.7, sigma = 0.4 and sqrt(0.21), and shapes 0.1
    to 3, the relative error is below 1e-14.

    method "series": the convergent series, as a check (descending for beta < 1, ascending for
    beta > 1, 1/(1 - i w) at beta = 1). A half is refused, with NaN and a warning, when its
    estimated rounding error 1e-15 sum_k |t_k| (1 + |log t_k| + |arg t_k|) exceeds tol relative to
    the half, when the series does not converge within 200000 terms, or when the sum breaks
    |G(w)| <= Gamma(1 + 1/beta). On scans of 400 arguments per shape for shapes 0.05 to 50 the
    accepted values agreed with the quadrature to about tol relative to each half (to within a
    factor 3 at tol = 1e-8); many arguments are refused, more for shapes near one."""
    _check_line(sigma, beta, lam1, lam2)
    if method not in ("quadrature", "series"):
        raise ValueError("method must be 'quadrature' or 'series'")
    tol = _check_tol(tol)
    n = len(_seq(t))
    arr = _d(_seq(t))
    re = (_ct.c_double * max(n, 1))()
    im = (_ct.c_double * max(n, 1))()
    _lib.wsgn_cf(arr, _i(n), _d1(mu), _d1(sigma), _d1(beta), _d1(lam1), _d1(lam2), _d1(tol),
                 _i(1 if method == "series" else 0), re, im)
    out = [complex(a, b) for a, b in zip(list(re)[:n], list(im)[:n])]
    if method == "series" and any(v != v for v in list(re)[:n]):
        import warnings as _w
        _w.warn("the series cancels beyond the requested precision at some arguments (NaN "
                "returned); use the default quadrature there")
    return out


def rsgn(n, mu=0.0, sigma=1.0, beta=2.0, lam1=0.5, lam2=0.5, wrap=False, rng=None):
    """Random draws by the uniform scale mixture: Z^beta ~ Gamma(1 + 1/beta), then uniform on
    [mu - Z sigma/lam2, mu + Z sigma/lam1]; reduced modulo 2 pi when wrap is true."""
    _check_line(sigma, beta, lam1, lam2)
    rng = rng or _random
    out = []
    for _ in range(n):
        u = rng.gammavariate(1.0 + 1.0 / beta, 1.0) ** (1.0 / beta)
        x = rng.uniform(mu - (sigma / lam2) * u, mu + (sigma / lam1) * u)
        out.append(x % _TWO_PI if wrap else x)
    return out


# ------------------------------------------------------------------ the wrapped distribution
def dwsgn(theta, mu=0.0, beta=2.0, lam=0.5, K=50, sigma=None, log=False):
    """Wrapped density (or log-density) at theta; parameters (mu, sigma, lam, beta), right rate
    lam/sigma, left rate (1-lam)/sigma; the four-parameter family is identifiable. sigma=None gives
    the reciprocal-rate subfamily sigma = sqrt(lam(1-lam)) of version 0.1.0, a three-parameter
    restriction. Images are summed directly on each side of mu and the remainder is added by the
    Euler-Maclaurin formula (exact incomplete-gamma integral, Bernoulli terms through B_12) as soon
    as it is accurate, tried after 8, 16, 32 and K images and after each further one, so the result
    does not depend on K. The tail is kept in log space, so very small shapes are covered (at beta
    = 0.005 and 0.001 the density differs from 1/(2 pi) by less than 1e-13 relative). Against
    40-digit image sums at six angles, lam = 0.3, sigma = 0.4 and sqrt(0.21), and shapes 0.1, 0.2,
    0.3, 0.5, 1 and 2, the relative error is below 1e-14."""
    s = _sigma(sigma, lam)
    _check(s, beta, lam, K)
    n = len(_seq(theta))
    arr = _d(_seq(theta))
    out = (_ct.c_double * max(n, 1))()
    _lib.wsgn_wdensity(arr, _i(n), _d1(mu), _d1(s), _d1(beta), _d1(lam), _i(K),
                       _i(1 if log else 0), out)
    return list(out)[:n]


def pwsgn(q, mu=0.0, beta=2.0, lam=0.5, K=50, sigma=None):
    """Wrapped distribution function, measured anticlockwise from angle 0: the probability of the
    arc [0, q]; 0 for q <= 0 and 1 for q >= 2 pi. Computed from the incomplete-gamma survival
    functions of the two sides summed over images, with the Euler-Maclaurin tail; the arc is
    described by its start and width, so a short arc keeps its relative accuracy (for an arc shorter
    than 0.01 the first image comes from the incomplete gamma function where that is accurate and
    from Gauss-Legendre quadrature of the density otherwise), and the incomplete gamma functions
    are evaluated from the logarithm of their argument, so (rate x)^beta never underflows next to
    the location. Small probabilities therefore keep their relative accuracy down to the underflow
    limit, on both sides of the location and for large shapes too (checked for shapes 0.3 to 50)."""
    s = _sigma(sigma, lam)
    _check(s, beta, lam, K)
    n = len(_seq(q))
    arr = _d(_seq(q))
    out = (_ct.c_double * max(n, 1))()
    _lib.wsgn_wcdf(arr, _i(n), _d1(mu), _d1(s), _d1(beta), _d1(lam), _i(K), out)
    return list(out)[:n]


def qwsgn(p, mu=0.0, beta=2.0, lam=0.5, K=50, sigma=None):
    """Wrapped quantile function on [0, 2 pi], by bisection on pwsgn with geometric midpoints
    while the bracket is wide, so that small probabilities get quantiles with relative accuracy
    (for instance qwsgn([1e-300]))."""
    s = _sigma(sigma, lam)
    _check(s, beta, lam, K)
    n = len(_seq(p))
    arr = _d(_seq(p))
    out = (_ct.c_double * max(n, 1))()
    _lib.wsgn_wquant(arr, _i(n), _d1(mu), _d1(s), _d1(beta), _d1(lam), _i(K), out)
    return list(out)[:n]


def rwsgn(n, mu=0.0, beta=2.0, lam=0.5, sigma=None, rng=None):
    """Draws from the wrapped distribution (angles in [0, 2 pi))."""
    s = _sigma(sigma, lam)
    _check(s, beta, lam)
    return rsgn(n, mu, s, beta, lam, 1.0 - lam, wrap=True, rng=rng)


def trig_moments(P=2, mu=0.0, beta=2.0, lam=0.5, sigma=None, tol=1e-12):
    try:
        pf = float(P)
    except (TypeError, ValueError):
        pf = float("nan")
    if not (_math.isfinite(pf) and pf == int(pf) and 1 <= pf <= 1e7):
        raise ValueError("'P' must be one whole number between 1 and 1e7")
    P = int(pf)
    """Trigonometric moments phi_p = A_p + i B_p, p = 1..P, of the wrapped distribution: the
    characteristic function at the integers, computed as in cfsgn, so their error is absolute, of
    order tol / p and in practice of order tol: about 1e-15 to 3e-15 for shapes up to about 3, a few times 1e-15 for shapes 5 to 20,
    and up to about 1e-14 for shapes up to 50, at the default tol; moments much smaller than that carry no relative
    accuracy. P is a whole number from 1 to 1e7; tol as in cfsgn."""
    s = _sigma(sigma, lam)
    _check(s, beta, lam)
    A = (_ct.c_double * max(P, 1))()
    B = (_ct.c_double * max(P, 1))()
    _lib.wsgn_trigmom(_i(P), _d1(mu), _d1(s), _d1(beta), _d1(lam), _d1(_check_tol(tol)), A, B)
    return {"A": list(A)[:P], "B": list(B)[:P]}


def trig_decay(beta=2.0, lam=0.5, sigma=None):
    """Constant C and exponent 1 + beta of the exact decay rate |phi_p| ~ C p^-(1 + beta), which
    holds unless lam = 0.5 and beta is an even integer (then C = 0)."""
    s = _sigma(sigma, lam)
    _check(s, beta, lam)
    out = (_ct.c_double * 2)()
    _lib.wsgn_decay(_d1(s), _d1(beta), _d1(lam), out)
    return {"C": out[0], "exponent": out[1]}


def nll(theta, mu, beta, lam, K=50, sigma=None):
    """Negative log-likelihood of a wrapped sample."""
    theta = [float(v) for v in theta]
    s = _sigma(sigma, lam)
    _check(s, beta, lam, K)
    arr = _d(theta)
    out = (_ct.c_double * 1)()
    _lib.wsgn_nll(arr, _i(len(theta)), _d1(mu), _d1(s), _d1(beta), _d1(lam), _i(K), out)
    return out[0]


# ------------------------------------------------------------------ the 0.1.0 fitters (subfamily)
def _nelder_mead(f, x0, maxit=500, tol=1e-8):
    n = len(x0)
    al, ga, rh, sg = 1.0, 2.0, 0.5, 0.5
    pts = [list(x0)]
    for i in range(n):
        p = list(x0)
        p[i] += (0.1 * p[i] if p[i] != 0 else 0.1)
        pts.append(p)
    vals = [f(p) for p in pts]
    for _ in range(maxit):
        order = sorted(range(n + 1), key=lambda k: vals[k])
        pts = [pts[k] for k in order]
        vals = [vals[k] for k in order]
        if abs(vals[-1] - vals[0]) < tol:
            break
        cen = [sum(p[i] for p in pts[:-1]) / n for i in range(n)]
        xr = [cen[i] + al * (cen[i] - pts[-1][i]) for i in range(n)]
        fr = f(xr)
        if fr < vals[0]:
            xe = [cen[i] + ga * (xr[i] - cen[i]) for i in range(n)]
            fe = f(xe)
            pts[-1], vals[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < vals[-2]:
            pts[-1], vals[-1] = xr, fr
        else:
            xc = [cen[i] + rh * (pts[-1][i] - cen[i]) for i in range(n)]
            fc = f(xc)
            if fc < vals[-1]:
                pts[-1], vals[-1] = xc, fc
            else:
                for k in range(1, n + 1):
                    pts[k] = [pts[0][i] + sg * (pts[k][i] - pts[0][i]) for i in range(n)]
                    vals[k] = f(pts[k])
    b = min(range(n + 1), key=lambda k: vals[k])
    return pts[b], vals[b]


def _solve3(a, b):
    """Solve the 3 x 3 system a x = b (b a list of columns) by Gauss-Jordan with pivoting."""
    m = [list(a[r]) + [col[r] for col in b] for r in range(3)]
    for c in range(3):
        piv = max(range(c, 3), key=lambda r: abs(m[r][c]))
        m[c], m[piv] = m[piv], m[c]
        d = m[c][c]
        m[c] = [v / d for v in m[c]]
        for r in range(3):
            if r != c:
                f = m[r][c]
                m[r] = [vr - f * vc for vr, vc in zip(m[r], m[c])]
    return [[m[r][3 + j] for r in range(3)] for j in range(len(b))]


def _md_sandwich_se(theta, est, mom):
    """Minimum-distance sandwich standard errors, as in the R front:
    (G'G)^{-1} G' V G (G'G)^{-1} / n, G by central differences, V the moment covariance."""
    n = len(theta)
    z = [[_math.cos(t), _math.sin(t), _math.cos(2 * t), _math.sin(2 * t)] for t in theta]
    zb = [sum(r[j] for r in z) / n for j in range(4)]
    V = [[sum((r[i] - zb[i]) * (r[j] - zb[j]) for r in z) / (n - 1) for j in range(4)]
         for i in range(4)]
    h = 1e-4
    G = []
    for j in range(3):
        e = [0.0, 0.0, 0.0]
        e[j] = h
        up = mom([est[k] + e[k] for k in range(3)])
        dn = mom([est[k] - e[k] for k in range(3)])
        G.append([(u - d) / (2 * h) for u, d in zip(up, dn)])
    GtG = [[sum(G[i][r] * G[j][r] for r in range(4)) for j in range(3)] for i in range(3)]
    cols = _solve3(GtG, [[G[i][r] for i in range(3)] for r in range(4)])
    Bm = [[cols[r][i] for r in range(4)] for i in range(3)]
    se = []
    for i in range(3):
        v = sum(Bm[i][a] * V[a][b] * Bm[i][b] for a in range(4) for b in range(4))
        se.append(_math.sqrt(v / n))
    return se


def fit_wsgn_trig(theta):
    """Minimum-distance trigonometric-moment fit of the three-parameter reciprocal-rate subfamily
    (sigma = sqrt(lam(1-lam))), kept from version 0.1.0; returns a dict with par (mu, beta, lam),
    the sandwich standard errors se, and the objective value. Useful as a diagnostic and for
    starting values; fit_wsgn fits the full four-parameter distribution."""
    n = len(theta)
    m = [sum(_math.cos(t) for t in theta) / n, sum(_math.sin(t) for t in theta) / n,
         sum(_math.cos(2 * t) for t in theta) / n, sum(_math.sin(2 * t) for t in theta) / n]
    mu0 = _math.atan2(m[1], m[0]) % _TWO_PI

    def mom(p3):
        r = trig_moments(2, p3[0], p3[1], p3[2])
        return [r["A"][0], r["B"][0], r["A"][1], r["B"][1]]

    def obj(p3):
        if p3[1] < 0.3 or p3[1] > 12 or p3[2] <= 0.02 or p3[2] >= 0.98:
            return 1e7
        return sum((a - b) ** 2 for a, b in zip(mom(p3), m))
    best = None
    for s in ([mu0, 2.0, 0.5], [mu0, 1.0, 0.3], [mu0, 4.0, 0.7]):
        p, v = _nelder_mead(obj, s)
        if best is None or v < best[1]:
            best = (p, v)
    se = _md_sandwich_se(theta, best[0], mom)
    names = ("mu", "beta", "lam")
    return {"par": dict(zip(names, best[0])), "se": dict(zip(names, se)), "value": best[1]}


def fit_wsgn_mle(theta, start=None, K=50, B=0, rng=None):
    """Maximum likelihood fit of the three-parameter reciprocal-rate subfamily, kept from version
    0.1.0; returns a dict with par and logLik, and, when B > 0, parametric-bootstrap standard
    errors se. fit_wsgn fits the full four-parameter distribution."""
    if start is None:
        s = fit_wsgn_trig(theta)["par"]
        start = [s["mu"], s["beta"], s["lam"]]

    def make_f(th):
        def f(p3):
            if p3[1] < 0.3 or p3[1] > 12 or p3[2] <= 0.02 or p3[2] >= 0.98:
                return 1e7
            v = nll(th, p3[0], p3[1], p3[2], K)
            return v if v == v and v not in (float("inf"), float("-inf")) else 1e7
        return f
    p, v = _nelder_mead(make_f(theta), list(start))
    names = ("mu", "beta", "lam")
    out = {"par": dict(zip(names, p)), "logLik": -v}
    if B > 0:
        ests = []
        for _ in range(B):
            x = rwsgn(len(theta), p[0], p[1], p[2], rng=rng)
            ests.append(_nelder_mead(make_f(x), list(p))[0])
        se = []
        for j in range(3):
            col = [e[j] for e in ests]
            mb = sum(col) / B
            se.append(_math.sqrt(sum((c - mb) ** 2 for c in col) / (B - 1)) if B > 1 else 0.0)
        out["se"] = dict(zip(names, se))
    return out


# ------------------------------------------------------------------ four-parameter likelihood
_MAX_CELLS = 10 ** 6
_ONE_SIDED_GAP = _math.log(1e12)


def _info_ok(info, vc, s4, sym):
    # finite, positive definite information (lam row and column dropped when symmetric)
    act = [0, 1, 3] if sym else [0, 1, 2, 3]
    I = [[info[i + 4 * j] for j in act] for i in act]
    if not all(_math.isfinite(v) for row in I for v in row):
        return False
    if not all(_math.isfinite(s4[i]) and s4[i] > 0 for i in act):
        return False
    m = len(act)
    L = [[0.0] * m for _ in range(m)]
    for i in range(m):
        for j in range(i + 1):
            t = I[i][j] - sum(L[i][k] * L[j][k] for k in range(j))
            if i == j:
                if not t > 0:
                    return False
                L[i][i] = _math.sqrt(t)
            else:
                L[i][j] = t / L[j][j]
    return True


def _fit_args(theta, beta_floor, resolution, origin, symmetric, window, K):
    theta = [float(v) for v in theta]
    if len(theta) < 3 or not all(_math.isfinite(v) for v in theta):
        raise ValueError("'theta' must be at least three finite angles (no NaN or inf)")
    if not _math.isfinite(origin):
        raise ValueError("'origin' must be one finite number")
    grouped = resolution is not None
    if grouped:
        if not (0 < resolution <= _math.pi):
            raise ValueError("'resolution' must be one positive number of radians, at most pi")
        J = round(_TWO_PI / resolution)
        if J > _MAX_CELLS:
            raise ValueError("'resolution' gives %d cells; at most 1e6 are allowed (resolution "
                             "at least 2 * pi / 1e6)" % J)
        if abs(J * resolution - _TWO_PI) > 1e-8:
            raise ValueError("'resolution' must divide the circle into a whole number of cells")
        if len(set(round((v - origin) / resolution) % J for v in theta)) < 2:
            raise ValueError("the data occupy a single cell: the grouped likelihood has no "
                             "interior maximum")
    elif len(set(v % _TWO_PI for v in theta)) < 3:
        raise ValueError("exact data need at least three distinct angles: with fewer, the "
                         "likelihood is unbounded at every shape")
    if beta_floor is None:
        beta_floor = 1e-3 if grouped else 0.1
    if not 0 < beta_floor < 20:
        raise ValueError("'beta_floor' must be one number in (0, 20)")
    _check_K(K)
    return dict(theta=theta, n=len(theta), grouped=1 if grouped else 0,
                res=float(resolution) if grouped else 0.0, origin=float(origin),
                beta0=float(beta_floor), sym=1 if symmetric else 0, window=float(window), K=int(K))


_OUT = ("mu", "sigma", "lam", "beta", "logLik", "nfev", "at_floor", "starts", "at_upper",
        "at_rate_bound", "log_a", "log_b")


_MAX_CORES = 256


def _cap_cores(ncores):
    if int(float(ncores)) > _MAX_CORES:
        import warnings as _w
        _w.warn("'ncores' reduced from %d to %d, the most threads used" % (int(float(ncores)), _MAX_CORES),
                stacklevel=3)
        return _MAX_CORES
    return int(float(ncores))


def _fit_batch(draws, a, B, par, ncores=1, lrates=None):
    """Refit each block of n draws (B blocks, block-major) by the C procedure; the parity harness
    passes the same draws to both front ends. lrates = (log a, log b) of the fit (from par if
    omitted)."""
    if lrates is None:
        lrates = (_math.log(par["lam"] / par["sigma"]), _math.log((1.0 - par["lam"]) / par["sigma"]))
    est = (_ct.c_double * (12 * B))()
    se = (_ct.c_double * 5)()
    _lib.wsgn_fit_batch(_d(draws), _i(a["n"]), _i(B), _i(a["grouped"]), _d1(a["res"]),
                        _d1(a["origin"]), _d1(a["beta0"]), _i(a["sym"]), _d1(a["window"]),
                        _i(a["K"]), _i(ncores), _d([par[k] for k in _NAMES] + list(lrates)), est,
                        se)
    rows = [dict(zip(_OUT[:5], est[12 * b:12 * b + 5])) for b in range(B)]
    return {"est": rows, "se": dict(zip(_NAMES, se[0:4])), "frac_beta_lt1": se[4]}


def fit_wsgn(theta, beta_floor=None, resolution=None, origin=0.0, symmetric=False, se="auto",
             B=200, window=0.4, K=50, ngrid=20000, ncores=1, rng=None):
    """Maximum likelihood fit of the four-parameter wrapped distribution (mu, sigma, lam, beta);
    every step (optimiser, starting values, location search, bootstrap refits) runs in C.

    Exact data: the estimate is the maximiser of the likelihood over beta_floor <= beta <= 20 and
    |log a|, |log b| <= R, R = max(30, 3 (1 - log beta_floor) / beta_floor) (99 at the default
    floor 0.1). Two searches are made and the best result is returned. The regular branch:
    Nelder-Mead from 27 starting values about the circular mean and 6 one-sided starts (one
    log-rate near the bound, the location at either end of the widest gap between observations),
    with the location searched over the stored observations within `window` of the fit when its
    shape is at most one. The spike branch, whose likelihood grows as the shape falls and whose
    profile can be ragged just above the floor: the same procedure, with a spike start on the
    most-tied observation, at the fixed shapes beta_floor times 1, 1.25 and 1.6 (the 1.25 and 1.6 shapes are skipped when
    they reach half the upper bound, i.e. when the floor is 6.25 or more), with the location
    over every observation, each followed by a free polish. In every location search each
    candidate observation is screened after a short optimisation of the other parameters (from the
    current fit and from a spike as wide as the gap to its nearest neighbour): all candidates when
    there are at most 200, otherwise the 15 best at the current fit, the 15 with the smallest gap
    to a neighbour, the 15 with the smallest gap to the further neighbour and every tied value; the
    12 best and every tie are refined further, and the best two fully. profile_wsgn uses the same
    search, so the profile at the floor equals the fit's floor search. at_floor, at_upper (upper
    bound 20) and at_rate_bound (the rate bound R: a degenerate fit, where lam may round to 0 or 1;
    log_rates are exact) report fits on the bounds; at_one_sided reports a one-sided fit that stops
    short of the rate bound (one rate effectively infinite, lam within about 1e-12 of 0 or 1, the
    log-rates differing by more than log(1e12)), for exact and grouped data alike, the true maximum on the boundary, kept as such. The search is not exhaustive: on the cases
    checked it matched or beat independent searches, but more distant branches cannot be excluded.
    Exact data need at least three distinct angles, all finite. Without standard errors a fit takes
    about 2.5 seconds for the 31 periwinkles, 5 to 10 seconds for 60 observations, 20 for 150 and
    40 for 300 on a current laptop.

    Rounded data: give `resolution` (radians); the grouped likelihood over cells of that width
    centred at origin + j*resolution (at most 1e6 cells; a finer resolution is refused) is
    maximised (default floor 0.001, upper bound 50, rate bound
    60), from 30 fixed-seed random starts; the data must occupy at least two cells.

    Standard errors: se="auto" gives the expected Fisher information for a continuous fit with
    beta > 1 on none of the bounds (in a simulation at n = 200 the 95 per cent Wald intervals
    covered 87 to 92 per cent), and otherwise no standard errors, with a warning saying why (on a
    bound they have no meaning; at a shape at most one or for grouped data the bootstrap is
    available on request). A one-sided fit (at_one_sided) counts as on a bound. Standard errors
    are never returned as NaN, except the lam entry of a symmetric fit, where lam is fixed: if the information matrix is not finite and positive definite,
    se="auto" gives none, with a warning. se="information" is refused for grouped data, a fit on
    a bound (the floor, the upper bound or the rate bound), a one-sided fit, and when the
    information matrix cannot be inverted.
    se="bootstrap" draws B samples from the fitted rates (rounded to the same cells when grouped)
    and refits each by the same procedure, as long as a fit each.
    Speed: ncores (a whole number, default 1) sets the number of threads when the library was built
    with OpenMP (HAS_OPENMP is True: the Linux wheels and a source build on Linux; the macOS wheels
    are serial). The independent trials inside a fit (the starting values, the candidate locations,
    the spike branches) then run in parallel, and with se="bootstrap" the refits do, each on one
    thread. The best trial is chosen afterwards in a fixed order, so the fit is identical to the
    last bit for every ncores. A likelihood evaluation with 128 or more observations (or cells) is
    itself split over the threads and summed in a fixed order. On a laptop with four performance
    and six efficiency cores, busy with other work, a fit of 1500 observations took 160 seconds
    on one thread, 86 on two, 57 on four and 55 on eight; the gain stops near four threads because
    the first pass of the location search, the Nelder-Mead runs and the final polish still run in
    sequence and the extra threads run on slower cores. A fit takes minutes from about a thousand
    observations on. ncores above 256 is reduced to 256 with a warning. Without OpenMP ncores is
    accepted and has no effect.
    Returns a dict with par, rates, log_rates, logLik, se, vcov, boot, frac_beta_lt1, se_method,
    beta_floor,
    at_floor, at_upper, at_rate_bound, at_one_sided, n, grouped, resolution, symmetric and nfev."""
    if se not in ("auto", "information", "bootstrap", "none"):
        raise ValueError("se must be 'auto', 'information', 'bootstrap' or 'none'")
    a = _fit_args(theta, beta_floor, resolution, origin, symmetric, window, K)

    def _whole(v, lo):
        try:
            f = float(v)
        except (TypeError, ValueError):
            return False
        return _math.isfinite(f) and f == int(f) and f >= lo
    if not (_whole(B, 2) and _whole(ngrid, 100) and _whole(ncores, 1)):
        raise ValueError("'B' (at least 2), 'ngrid' (at least 100) and 'ncores' must be whole numbers")
    if float(ngrid) > 1e8:
        raise ValueError("'ngrid' must be at most 1e8")
    ncores = _cap_cores(ncores)
    if not (_math.isfinite(float(window)) and float(window) >= 0):
        raise ValueError("'window' must be one non-negative number")
    out8 = (_ct.c_double * 12)()
    _lib.wsgn_fit(_d(a["theta"]), _i(a["n"]), _i(a["grouped"]), _d1(a["res"]), _d1(a["origin"]),
                  _d1(a["beta0"]), _i(a["sym"]), _d1(0.0), _d1(a["window"]), _i(a["K"]), _i(ncores),
                  out8)
    r = dict(zip(_OUT, list(out8)))
    par = {k: r[k] for k in _NAMES}
    # one-sided: min(lam, 1 - lam) below about 1e-12 (one rate effectively infinite), judged from
    # the exact log-rates, which stay valid when lam rounds to 0 or 1
    one_sided = abs(r["log_a"] - r["log_b"]) > _ONE_SIDED_GAP
    on_bound = bool(r["at_floor"]) or bool(r["at_upper"]) or bool(r["at_rate_bound"])
    bound = on_bound or one_sided
    method = se
    if se == "auto":
        method = ("information" if (not a["grouped"] and par["beta"] > 1 and not bound)
                  else "none")
        if method == "none":
            import warnings as _w
            why = ("the fit is one-sided (lam is numerically 0 or 1, one rate effectively infinite), "
                   "on the boundary of the parameter space, where they have no meaning"
                   if (one_sided and not on_bound) else
                   "the fit is on a bound of the parameter space, where they have no meaning" if bound
                   else "the information is not used for grouped data" if a["grouped"]
                   else "the fitted shape is at most one, where the normal approximation is not "
                        "proved")
            _w.warn("no standard errors: " + why +
                    ("" if bound else "; se='bootstrap' gives descriptive bootstrap standard "
                     "errors"))
    if method == "information" and (a["grouped"] or bound):
        raise ValueError("information standard errors are not available for grouped data, for a fit "
                         "on a bound or for a one-sided fit (lam numerically 0 or 1); use "
                         "se='bootstrap'")
    res = {"par": par, "rates": {"a": _math.exp(r["log_a"]), "b": _math.exp(r["log_b"])},
           "log_rates": {"log_a": r["log_a"], "log_b": r["log_b"]},
           "logLik": r["logLik"], "se": None, "vcov": None, "boot": None, "frac_beta_lt1": None,
           "se_method": method, "beta_floor": a["beta0"], "at_floor": bool(r["at_floor"]),
           "at_upper": bool(r["at_upper"]), "at_rate_bound": bool(r["at_rate_bound"]),
           "at_one_sided": bool(one_sided),
           "n": a["n"], "grouped": bool(a["grouped"]), "resolution": resolution,
           "symmetric": bool(symmetric), "nfev": r["nfev"]}
    if method == "information":
        info = (_ct.c_double * 16)()
        vc = (_ct.c_double * 16)()
        s4 = (_ct.c_double * 4)()
        _lib.wsgn_info(_d([par[k] for k in _NAMES]), _i(a["sym"]), _i(a["K"]), _i(ngrid),
                       _i(a["n"]), info, vc, s4)
        if _info_ok(info, vc, s4, a["sym"]):
            res["se"] = dict(zip(_NAMES, list(s4)))
            res["vcov"] = [[vc[i + 4 * j] for j in range(4)] for i in range(4)]
        elif se == "information":
            raise ValueError("the information matrix is not finite and positive definite at this "
                             "fit, so there are no information standard errors; use "
                             "se='bootstrap'")
        else:
            import warnings as _w
            _w.warn("no standard errors: the information matrix is not finite and positive "
                    "definite at this fit; se='bootstrap' gives descriptive bootstrap standard "
                    "errors")
            res["se_method"] = "none"
    elif method == "bootstrap":
        # drawn from the rates, which stay valid when lam rounds to 0 or 1 on a degenerate fit
        draws = rsgn(a["n"] * int(B), par["mu"], 1.0, par["beta"], res["rates"]["a"],
                     res["rates"]["b"], wrap=True, rng=rng)
        bt = _fit_batch(draws, a, int(B), par, ncores, (r["log_a"], r["log_b"]))
        res["se"] = bt["se"]
        res["boot"] = bt["est"]
        res["frac_beta_lt1"] = bt["frac_beta_lt1"]
    return res


def profile_wsgn(theta, beta, resolution=None, origin=0.0, symmetric=False, window=_math.pi,
                 K=50, beta_floor=None, ncores=1):
    """Profile log-likelihood M(beta): the supremum over (mu, sigma, lam) at each fixed shape, for
    exact data or (with resolution) the grouped likelihood. For exact data it tends to infinity as
    beta tends to zero, at the leading rate n log(n/(n - m))/beta, m the largest number of tied
    observations. At shapes up to one the location is searched over the stored observations within
    `window` of the best free fit, as described for fit_wsgn (every observation screened when there
    are at most 200); the default window = pi covers all of them. The rate bound is that of fit_wsgn
    with the same beta_floor, max(30, 3 (1 - log b) / b) with b the smaller of beta_floor and the
    shape (at or above the floor, the bound is computed from beta_floor, not from each profiled
    shape), so the profile at the floor equals the floor search of the fit. ncores sets the
    threads, as for fit_wsgn (the shapes and the trials inside each fit share them); the profile
    is identical for every ncores, and ncores has no effect without OpenMP. Returns a list of dicts
    with keys beta, logLik, mu, sigma, lam."""
    try:
        ok = _math.isfinite(float(ncores)) and float(ncores) == int(ncores) and float(ncores) >= 1
    except (TypeError, ValueError):
        ok = False
    if not ok:
        raise ValueError("'ncores' must be a whole number, at least 1")
    ncores = _cap_cores(ncores)
    betas = [float(b) for b in _seq(beta)]
    if any(not (b > 0 and _math.isfinite(b)) for b in betas):
        raise ValueError("'beta' must be positive and finite")
    a = _fit_args(theta, beta_floor, resolution, origin, symmetric, window, K)
    nb = len(betas)
    out = (_ct.c_double * (12 * nb))()
    _lib.wsgn_profile(_d(a["theta"]), _i(a["n"]), _i(a["grouped"]), _d1(a["res"]),
                      _d1(a["origin"]), _d1(a["beta0"]), _d(betas), _i(nb), _i(a["sym"]),
                      _d1(a["window"]), _i(a["K"]), _i(ncores), out)
    return [{"beta": betas[i], "logLik": out[12 * i + 4], "mu": out[12 * i],
             "sigma": out[12 * i + 1], "lam": out[12 * i + 2]} for i in range(nb)]


def grouped_loglik(theta, resolution, mu, sigma, lam, beta, origin=0.0, K=50):
    """Grouped log-likelihood of rounded circular data: theta is assigned to the cells of width
    `resolution` centred at origin + j*resolution, and the result is sum_j N_j log P(Theta in C_j),
    with the cell probabilities computed exactly as in pwsgn. It is at most zero and bounded for
    every parameter value; fit_wsgn maximises it when resolution is given. Angles must be finite."""
    if resolution is None:
        raise ValueError("'resolution' must be given, as one number of radians")
    resolution = float(resolution)
    theta = [float(v) for v in theta]
    if not all(_math.isfinite(v) for v in theta) or not _math.isfinite(origin):
        raise ValueError("'theta' and 'origin' must be finite (no NaN or inf)")
    if not (0 < resolution <= _math.pi):
        raise ValueError("'resolution' must divide the circle into a whole number of cells")
    if round(_TWO_PI / resolution) > _MAX_CELLS:
        raise ValueError("'resolution' gives %d cells; at most 1e6 are allowed (resolution at "
                         "least 2 * pi / 1e6)" % round(_TWO_PI / resolution))
    if abs(round(_TWO_PI / resolution) * resolution - _TWO_PI) > 1e-8:
        raise ValueError("'resolution' must divide the circle into a whole number of cells")
    _check(sigma, beta, lam, K)
    a = dict(theta=theta, n=len(theta), res=float(resolution), origin=float(origin), K=int(K))
    out = (_ct.c_double * 1)()
    _lib.wsgn_grouped_ll(_d(a["theta"]), _i(a["n"]), _d1(a["res"]), _d1(a["origin"]), _d1(mu),
                         _d1(sigma), _d1(beta), _d1(lam), _i(K), out)
    return out[0]
