"""Parity harness, Python front: prints every value of every entry point of the shared C back-end as
"label %.15e", in the same order as parity_R.R. Run through parity_run.sh."""
import math
import wsgnstat as w
from wsgnstat import _core as C


def out(label, v):
    for i, x in enumerate(v):
        print("%s[%d] %.15e" % (label, i, x))


def cplx(label, z):
    out(label + ".re", [c.real for c in z])
    out(label + ".im", [c.imag for c in z])


x = [-3, -0.7, 0, 0.4, 1.9, 5]
th = [0, 0.3, 1.1, 2.5, 3.9, 5.2, 6.2]
tt = [0, 0.5, 1, 2, 3, 7.5, 20]
shapes = [0.05, 0.1, 0.3, 0.7, 1, 1.5, 2, 3.5, 8]
import warnings
warnings.simplefilter("ignore")
for be in shapes:
    lab = "b%.2f" % be
    out("dsgn." + lab, w.dsgn(x, 0.2, 0.8, be, 0.3, 0.7))
    out("psgn." + lab, w.psgn(x, 0.2, 0.8, be, 0.3, 0.7))
    cplx("cfsgn." + lab, w.cfsgn(tt, 0.2, 0.8, be, 0.3, 0.7))
    cplx("cfsgn.tol8." + lab, w.cfsgn(tt, 0.2, 0.8, be, 0.3, 0.7, tol=1e-8))
    cplx("cfser." + lab, w.cfsgn(tt, 0.2, 0.8, be, 0.3, 0.7, method="series"))
    tm = w.trig_moments(6, 1.3, be, 0.35, sigma=0.5)
    out("trigA." + lab, tm["A"]); out("trigB." + lab, tm["B"])
    tm0 = w.trig_moments(3, 1.3, be, 0.35)
    out("trig0A." + lab, tm0["A"]); out("trig0B." + lab, tm0["B"])
    dc = w.trig_decay(be, 0.35, 0.5)
    out("decay." + lab, [dc["C"], dc["exponent"]])
    out("dwsgn." + lab, w.dwsgn(th, 1.3, be, 0.35, sigma=0.5))
    out("dwsgnlog.K3." + lab, w.dwsgn(th, 1.3, be, 0.35, K=3, sigma=0.5, log=True))
    out("dwsgn0." + lab, w.dwsgn(th, 1.3, be, 0.35))
    out("pwsgn." + lab, w.pwsgn(th, 1.3, be, 0.35, sigma=0.5))
    out("qwsgn." + lab, w.qwsgn([0.01, 0.25, 0.5, 0.9], 1.3, be, 0.35, sigma=0.5))
    out("nll." + lab, [w.nll(th, 1.3, be, 0.35, sigma=0.5)])
    out("gll." + lab, [w.grouped_loglik(th, math.pi / 18, 1.3, 0.5, 0.35, be)])


def vdc(n):
    res = []
    for i in range(1, n + 1):
        s, f, k = 0.0, 0.5, i
        while k > 0:
            s += f * (k % 2); k //= 2; f /= 2
        res.append(s)
    return res


NAMES = ("mu", "sigma", "lam", "beta")
s1 = w.qwsgn(vdc(120), 1, 1.6, 0.3, sigma=0.4)
s2 = w.qwsgn(vdc(90), 2, 0.6, 0.4, sigma=0.25)
h = math.pi / 18
s3 = [(round(v / h) * h) % (2 * math.pi) for v in w.qwsgn(vdc(100), 3, 1.2, 0.45, sigma=0.3)]


def fo(f):
    v = [f["par"][k] for k in NAMES] + [f["logLik"], f["nfev"], float(f["at_floor"])]
    if f["se"] is not None:
        v += [f["se"][k] for k in NAMES]
    return v


f1 = w.fit_wsgn(s1)
out("fit.s1.info", fo(f1))
out("fit.s1.vcov", [f1["vcov"][i][j] for j in range(4) for i in range(4)])
out("fit.s1.sym", fo(w.fit_wsgn(s1, symmetric=True)))
out("fit.s2", fo(w.fit_wsgn(s2, se="none")))
out("fit.s2.floor", fo(w.fit_wsgn(s2, beta_floor=1.2, se="none")))
out("fit.s3.grouped", fo(w.fit_wsgn(s3, resolution=h, se="none")))
pr = w.profile_wsgn(s2, [0.1, 0.5, 1, 2])
out("profile.s2", [r[k] for k in ("beta", "logLik", "mu", "sigma", "lam") for r in pr])
pr3 = w.profile_wsgn(s3, [0.3, 1], resolution=h)
out("profile.s3", [r[k] for k in ("beta", "logLik", "mu", "sigma", "lam") for r in pr3])
a2 = C._fit_args(s2, None, None, 0.0, False, 0.4, 50)
f2 = w.fit_wsgn(s2, se="none")
p2 = f2["par"]
d2 = w.qwsgn(vdc(90 * 3), p2["mu"], p2["beta"], p2["lam"], sigma=p2["sigma"])
b2 = C._fit_batch(d2, a2, 3, p2)
out("batch.s2.est", [r[k] for r in b2["est"] for k in ("mu", "sigma", "lam", "beta", "logLik")])
out("batch.s2.se", [b2["se"][k] for k in NAMES] + [b2["frac_beta_lt1"]])
a3 = C._fit_args(s3, None, h, 0.0, False, 0.4, 50)
f3 = w.fit_wsgn(s3, resolution=h, se="none")
p3 = f3["par"]
d3 = w.qwsgn(vdc(100 * 2), p3["mu"], p3["beta"], p3["lam"], sigma=p3["sigma"])
b3 = C._fit_batch(d3, a3, 2, p3)
out("batch.s3.est", [r[k] for r in b3["est"] for k in ("mu", "sigma", "lam", "beta", "logLik")])
out("batch.s3.se", [b3["se"][k] for k in NAMES] + [b3["frac_beta_lt1"]])
