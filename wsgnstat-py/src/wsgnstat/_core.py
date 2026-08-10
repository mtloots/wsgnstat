"""wsgnstat: Python (ctypes) binding to the shared skew-generalised-normal C back-end. Pure
standard library. The C source is compiled on first import; the same source backs the R package
'wsgnstat'. The characteristic function is evaluated through absolutely convergent integrals (the
Fox H route), so densities, trigonometric moments and both estimators are available for every shape
parameter beta > 0, not only where the Fox-Wright series converges."""
import ctypes as _ct, math as _math, os as _os, random as _random

_here = _os.path.dirname(_os.path.abspath(__file__))
_ext = "dylib" if _os.uname().sysname == "Darwin" else "so"
_libpath = _os.path.join(_here, "libwsgnstat." + _ext)
if not _os.path.exists(_libpath):
    import subprocess as _sp
    _sp.check_call(["cc", "-O2", "-fPIC", "-shared", "-I", _here,
                    _os.path.join(_here, "wsgn.c"), "-lm", "-o", _libpath])
_lib = _ct.CDLL(_libpath)
_dp, _ip = _ct.POINTER(_ct.c_double), _ct.POINTER(_ct.c_int)
_lib.wsgn_density.argtypes  = [_dp, _ip, _dp, _dp, _dp, _dp, _dp, _dp]
_lib.wsgn_wdensity.argtypes = [_dp, _ip, _dp, _dp, _dp, _ip, _dp]
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

def dwsgn(theta, mu=0.0, beta=2.0, lam=0.5, K=50):
    """Wrapped density on [0, 2*pi), identified parameterisation (mu = location, beta, lam)."""
    arr = _d(theta); out = (_ct.c_double * len(arr))()
    _lib.wsgn_wdensity(arr, _i(len(arr)), _d1(mu), _d1(beta), _d1(lam), _i(K), out)
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


def fit_wsgn_trig(theta):
    """Minimum-distance trigonometric-moment fit; returns dict with par (mu, beta, lam) and value.
    The estimating equations match the R front; the simplex optimiser differs only in book-keeping,
    so estimates agree to optimiser tolerance rather than byte-identically."""
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
    return {"par": {"mu": best[0][0], "beta": best[0][1], "lam": best[0][2]}, "value": best[1]}


def fit_wsgn_mle(theta, start=None, K=50):
    """Maximum likelihood fit of the wrapped identified family; returns dict with par and logLik."""
    if start is None:
        s = fit_wsgn_trig(theta)["par"]
        start = [s["mu"], s["beta"], s["lam"]]
    def f(p3):
        if p3[1] < 0.3 or p3[1] > 12 or p3[2] <= 0.02 or p3[2] >= 0.98:
            return 1e7
        v = nll(theta, p3[0], p3[1], p3[2], K)
        return v if v == v and v not in (float("inf"), float("-inf")) else 1e7
    p, v = _nelder_mead(f, start)
    return {"par": {"mu": p[0], "beta": p[1], "lam": p[2]}, "logLik": -v}
