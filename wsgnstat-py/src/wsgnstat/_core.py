"""wsgnstat: Python (ctypes) binding to the shared skew-generalised-normal C back-end. Pure
standard library. The C source is compiled on first import; the same source backs the R package
'wsgnstat'. The characteristic function is evaluated through absolutely convergent integrals (the
Fox H route), so densities, trigonometric moments and both estimators are available for every shape
parameter beta > 0, not only where the Fox-Wright series converges."""
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
    try:
        _sp.check_call([cc, "-O2", "-fPIC", "-shared", "-I", _here] + srcs + ["-lm", "-o", out])
    except (OSError, _sp.CalledProcessError) as exc:
        raise ImportError(
            "wsgnstat could not build its C back-end from the bundled sources (%s). In an "
            "installed copy this should never happen, because pip builds the library at install "
            "time. In a source checkout a C compiler is required; set CC if it is not 'cc'." % exc)
    return out


_libpath = _find_built()
if _libpath is None:
    _cached = _os.path.join(_cache_dir(), "libwsgnstat." + _ext())
    _srcs = [f for f in _glob.glob(_os.path.join(_here, "*.c"))]
    _stale = (not _os.path.exists(_cached) or
              any(_os.path.getmtime(f) > _os.path.getmtime(_cached) for f in _srcs))
    _libpath = _compile_from_source() if _stale else _cached
if _libpath is None:
    raise ImportError("wsgnstat: no compiled back-end and no C sources to build one from.")

_lib = _ct.CDLL(_libpath)
_dp, _ip = _ct.POINTER(_ct.c_double), _ct.POINTER(_ct.c_int)
_lib.wsgn_density.argtypes  = [_dp, _ip, _dp, _dp, _dp, _dp, _dp, _dp]
_lib.wsgn_pdist.argtypes    = [_dp, _ip, _dp, _dp, _dp, _dp, _dp, _dp]
_lib.wsgn_wdensity.argtypes = [_dp, _ip, _dp, _dp, _dp, _ip, _dp]
_lib.wsgn_wcdf.argtypes     = [_dp, _ip, _dp, _dp, _dp, _ip, _dp]
_lib.wsgn_wquant.argtypes   = [_dp, _ip, _dp, _dp, _dp, _ip, _dp]
_lib.wsgn_cf.argtypes       = [_dp, _ip, _dp, _dp, _dp, _dp, _dp, _dp, _dp]
_lib.wsgn_trigmom.argtypes  = [_ip, _dp, _dp, _dp, _dp, _dp]
_lib.wsgn_nll.argtypes      = [_dp, _ip, _dp, _dp, _dp, _ip, _dp]


def _d(seq):
    seq = [float(v) for v in seq]
    a = (_ct.c_double * max(len(seq), 1))()
    for i, v in enumerate(seq):
        a[i] = v
    return a
def _i(v): a = (_ct.c_int * 1)(); a[0] = int(v); return a
def _d1(v): a = (_ct.c_double * 1)(); a[0] = float(v); return a


def dsgn(x, mu=0.0, sigma=1.0, beta=2.0, lam1=0.5, lam2=0.5):
    """Density of the skew generalised normal distribution on the line at x."""
    arr = _d(x); out = (_ct.c_double * len(arr))()
    _lib.wsgn_density(arr, _i(len(arr)), _d1(mu), _d1(sigma), _d1(beta), _d1(lam1), _d1(lam2), out)
    return list(out)

def psgn(q, mu=0.0, sigma=1.0, beta=2.0, lam1=0.5, lam2=0.5):
    """Distribution function of the line law at q, in closed form (regularised incomplete gamma)."""
    arr = _d(q); out = (_ct.c_double * len(arr))()
    _lib.wsgn_pdist(arr, _i(len(arr)), _d1(mu), _d1(sigma), _d1(beta), _d1(lam1), _d1(lam2), out)
    return list(out)

def dwsgn(theta, mu=0.0, beta=2.0, lam=0.5, K=50):
    """Wrapped density on [0, 2*pi), identified parameterisation (mu = location, beta, lam)."""
    arr = _d(theta); out = (_ct.c_double * len(arr))()
    _lib.wsgn_wdensity(arr, _i(len(arr)), _d1(mu), _d1(beta), _d1(lam), _i(K), out)
    return list(out)

def pwsgn(q, mu=0.0, beta=2.0, lam=0.5, K=50):
    """Wrapped distribution function, measured anticlockwise from angle 0: the probability of the
    arc [0, q]. For q >= 2*pi it returns the mass the wrapping sum truncated at K captures."""
    arr = _d(q); out = (_ct.c_double * len(arr))()
    _lib.wsgn_wcdf(arr, _i(len(arr)), _d1(mu), _d1(beta), _d1(lam), _i(K), out)
    return list(out)

def qwsgn(p, mu=0.0, beta=2.0, lam=0.5, K=50):
    """Wrapped quantile function on [0, 2*pi], by bisection on pwsgn. Warns for p above the mass
    the truncated wrapping sum captures."""
    tot = pwsgn([2.0 * _math.pi], mu, beta, lam, K)[0]
    if any(tot < float(v) < 1.0 for v in p):
        import warnings as _w
        _w.warn("the wrapping sum truncated at K = %d captures mass %g; raise K for probabilities "
                "above it" % (K, tot))
    arr = _d(p); out = (_ct.c_double * len(arr))()
    _lib.wsgn_wquant(arr, _i(len(arr)), _d1(mu), _d1(beta), _d1(lam), _i(K), out)
    return list(out)

def cfsgn(t, mu=0.0, sigma=1.0, beta=2.0, lam1=0.5, lam2=0.5):
    """Characteristic function of the line law at t, as a list of complex numbers."""
    arr = _d(t); re = (_ct.c_double * len(arr))(); im = (_ct.c_double * len(arr))()
    _lib.wsgn_cf(arr, _i(len(arr)), _d1(mu), _d1(sigma), _d1(beta), _d1(lam1), _d1(lam2), re, im)
    return [complex(a, b) for a, b in zip(re, im)]

def trig_moments(P=2, mu=0.0, beta=2.0, lam=0.5):
    """Population trigonometric moments (A_p, B_p), p = 1..P, of the wrapped family."""
    A = (_ct.c_double * P)(); B = (_ct.c_double * P)()
    _lib.wsgn_trigmom(_i(P), _d1(mu), _d1(beta), _d1(lam), A, B)
    return {"A": list(A), "B": list(B)}

def nll(theta, mu, beta, lam, K=50):
    """Negative log-likelihood of a wrapped sample under the identified parameterisation."""
    arr = _d(theta); out = (_ct.c_double * 1)()
    _lib.wsgn_nll(arr, _i(len(arr)), _d1(mu), _d1(beta), _d1(lam), _i(K), out)
    return out[0]

def rsgn(n, mu=0.0, sigma=1.0, beta=2.0, lam1=0.5, lam2=0.5, wrap=False, rng=None):
    """Random draws via the stochastic representation (U^beta gamma-distributed)."""
    rng = rng or _random
    out = []
    for _ in range(n):
        u = rng.gammavariate(1.0 + 1.0/beta, 1.0) ** (1.0/beta)
        x = rng.uniform(mu - (sigma/lam2)*u, mu + (sigma/lam1)*u)
        out.append(x % (2.0*_math.pi) if wrap else x)
    return out


def rwsgn(n, mu=0.0, beta=2.0, lam=0.5, rng=None):
    """Draws from the wrapped law in the identified parameterisation (angles in [0, 2*pi))."""
    return rsgn(n, mu, _math.sqrt(lam * (1.0 - lam)), beta, lam, 1.0 - lam, wrap=True, rng=rng)


def _nelder_mead(f, x0, maxit=500, tol=1e-8):
    n = len(x0); al, ga, rh, sg = 1.0, 2.0, 0.5, 0.5
    pts = [list(x0)]
    for i in range(n):
        p = list(x0); p[i] += (0.1 * p[i] if p[i] != 0 else 0.1); pts.append(p)
    vals = [f(p) for p in pts]
    for _ in range(maxit):
        order = sorted(range(n + 1), key=lambda k: vals[k])
        pts = [pts[k] for k in order]; vals = [vals[k] for k in order]
        if abs(vals[-1] - vals[0]) < tol:
            break
        cen = [sum(p[i] for p in pts[:-1]) / n for i in range(n)]
        xr = [cen[i] + al * (cen[i] - pts[-1][i]) for i in range(n)]; fr = f(xr)
        if fr < vals[0]:
            xe = [cen[i] + ga * (xr[i] - cen[i]) for i in range(n)]; fe = f(xe)
            pts[-1], vals[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < vals[-2]:
            pts[-1], vals[-1] = xr, fr
        else:
            xc = [cen[i] + rh * (pts[-1][i] - cen[i]) for i in range(n)]; fc = f(xc)
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
    G = []                                    # G[j] is column j (4 values)
    for j in range(3):
        e = [0.0, 0.0, 0.0]; e[j] = h
        up = mom([est[k] + e[k] for k in range(3)])
        dn = mom([est[k] - e[k] for k in range(3)])
        G.append([(u - d) / (2 * h) for u, d in zip(up, dn)])
    GtG = [[sum(G[i][r] * G[j][r] for r in range(4)) for j in range(3)] for i in range(3)]
    # Bm = (G'G)^{-1} G' is 3 x 4; columns of G' are the rows of G, i.e. G[.][r]
    cols = _solve3(GtG, [[G[i][r] for i in range(3)] for r in range(4)])   # 4 columns of Bm
    Bm = [[cols[r][i] for r in range(4)] for i in range(3)]
    se = []
    for i in range(3):
        v = sum(Bm[i][a] * V[a][b] * Bm[i][b] for a in range(4) for b in range(4))
        se.append(_math.sqrt(v / n))
    return se


def fit_wsgn_trig(theta):
    """Minimum-distance trigonometric-moment fit; returns a dict with par (mu, beta, lam), the
    sandwich standard errors se, and the objective value. The estimating equations match the R
    front; the simplex optimiser differs only in book-keeping, so estimates agree with R to
    optimiser tolerance rather than byte-identically."""
    n = len(theta)
    m = [sum(_math.cos(t) for t in theta)/n, sum(_math.sin(t) for t in theta)/n,
         sum(_math.cos(2*t) for t in theta)/n, sum(_math.sin(2*t) for t in theta)/n]
    mu0 = _math.atan2(m[1], m[0]) % (2.0*_math.pi)
    def mom(p3):
        r = trig_moments(2, p3[0], p3[1], p3[2])
        return [r["A"][0], r["B"][0], r["A"][1], r["B"][1]]
    def obj(p3):
        if p3[1] < 0.3 or p3[1] > 12 or p3[2] <= 0.02 or p3[2] >= 0.98:
            return 1e7
        return sum((a - b)**2 for a, b in zip(mom(p3), m))
    best = None
    for s in ([mu0, 2.0, 0.5], [mu0, 1.0, 0.3], [mu0, 4.0, 0.7]):
        p, v = _nelder_mead(obj, s)
        if best is None or v < best[1]:
            best = (p, v)
    se = _md_sandwich_se(theta, best[0], mom)
    names = ("mu", "beta", "lam")
    return {"par": dict(zip(names, best[0])), "se": dict(zip(names, se)), "value": best[1]}


def fit_wsgn_mle(theta, start=None, K=50, B=0, rng=None):
    """Maximum likelihood fit of the wrapped identified family; returns a dict with par and logLik,
    and, when B > 0, parametric-bootstrap standard errors se from B refits to samples drawn at the
    estimate (the observed information is unstable near the cusp at beta = 1)."""
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
