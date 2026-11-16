"""Tests of the Python front. They defend the same correctness claims as the R tests: the
characteristic function against 120-digit series values for shapes 0.1 to 3 and equal to one at
zero, closed forms at shapes one and two, the wrapped density against 40-digit image sums and
integrating to one, the distribution and quantile functions, the mixture representation by
simulation, the decay constant, and the four-parameter fits (shape floor, location search, grouped
likelihood, information and bootstrap standard errors)."""
import csv
import math
import os
import random

import pytest

import wsgnstat as w

TWO_PI = 2.0 * math.pi
HERE = os.path.dirname(os.path.abspath(__file__))


def simpson(f, a, b, n=20000):
    """Composite Simpson on [a, b] with n (even) panels, f vectorised over a list."""
    h = (b - a) / n
    xs = [a + h * i for i in range(n + 1)]
    ys = f(xs)
    s = ys[0] + ys[-1] + 4 * sum(ys[1:-1:2]) + 2 * sum(ys[2:-1:2])
    return s * h / 3.0


def ks_stat(u):
    u = sorted(u)
    n = len(u)
    return max(max((i + 1) / n - v, v - i / n) for i, v in enumerate(u))


def rows(name):
    with open(os.path.join(HERE, name)) as fh:
        return list(csv.DictReader(fh))


# ------------------------------------------------------------------ characteristic function
def test_cf_is_one_at_zero_for_every_shape():
    for be in (0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1, 1.5, 2, 3, 10):
        assert abs(w.cfsgn([0.0], 0, 0.4, be, 0.3, 0.7)[0] - 1) < 1e-14


def test_cf_matches_series_reference():
    worst = 0.0
    for r in rows("ref_cf.csv"):
        s = 0.4 if r["case"] == "1" else math.sqrt(0.21)
        z = w.cfsgn([float(r["t"])], 0, s, float(r["beta"]), 0.3, 0.7)[0]
        zr = complex(float(r["re"]), float(r["im"]))
        worst = max(worst, abs(z - zr) / abs(zr))
    assert worst < 1e-12


def test_series_method_agrees_and_refuses_when_cancelling():
    import warnings
    for be in (0.1, 0.3, 0.6, 1, 2, 3):
        t = [0.5, 1, 2, 3]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            a = w.cfsgn(t, 0.2, 0.4, be, 0.3, 0.7, method="series")
        b = w.cfsgn(t, 0.2, 0.4, be, 0.3, 0.7)
        ok = [(x, y) for x, y in zip(a, b) if x.real == x.real]
        assert len(ok) >= 2
        assert max(abs(x - y) for x, y in ok) < 1e-12
    with pytest.warns(UserWarning, match="cancels"):
        z = w.cfsgn([20.0], 0, 0.4, 1.2, 0.3, 0.7, method="series")
    assert z[0].real != z[0].real
    # a sum stopped inside a cancellation once returned about 5e215 here
    with pytest.warns(UserWarning, match="cancels"):
        z = w.cfsgn([5.0], 0, 0.6, 1.5, 0.8, 0.2, method="series")
    assert z[0].real != z[0].real
    # over a grid, the series either agrees with the quadrature to 1e-8 or is refused
    import itertools
    import warnings
    worst, kept, total = 0.0, 0, 0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for be, l, s in itertools.product((0.1, 0.3, 0.7, 1, 1.1, 1.5, 2, 3, 5),
                                          (0.1, 0.5, 0.8, 0.95), (0.05, 0.6, 2)):
            tt = [0.5, 2, 5, 30]
            zs = w.cfsgn(tt, 0.3, s, be, l, 1 - l, method="series")
            zq = w.cfsgn(tt, 0.3, s, be, l, 1 - l)
            for a_, b_ in zip(zs, zq):
                total += 1
                if a_.real == a_.real:
                    kept += 1
                    worst = max(worst, abs(a_ - b_))
    assert worst < 1e-8 and kept > 0.5 * total
    with pytest.raises(ValueError):
        w.cfsgn([1.0], tol=-1)


def test_cf_closed_forms():
    a, b = 0.3 / 0.4, 0.7 / 0.4
    for t in (0.3, 1, 2.5, 7, 40):
        z = w.cfsgn([t], 0.5, 0.4, 1, 0.3, 0.7)[0]
        zc = complex(math.cos(0.5 * t), math.sin(0.5 * t)) * a * b / ((a - 1j * t) * (b + 1j * t))
        assert abs(z - zc) < 1e-14
        assert abs(w.cfsgn([t], 0, 0.4, 2, 0.5, 0.5)[0] - math.exp(-0.16 * t * t)) < 1e-14


def test_trig_moments_anchor_and_decay():
    tm = w.trig_moments(5, 0, 2, 0.5, sigma=0.6)
    assert max(abs(a - math.exp(-0.36 * p * p)) for p, a in zip(range(1, 6), tm["A"])) < 1e-14
    assert w.trig_moments(2, 1, 1.5, 0.3) == w.trig_moments(2, 1, 1.5, 0.3, sigma=math.sqrt(0.21))
    for be in (0.8, 1.5, 3):
        dc = w.trig_decay(be, 0.3, 0.4)
        assert dc["exponent"] == 1 + be
        tm = w.trig_moments(500, 0, be, 0.3, sigma=0.4)
        r = [abs(complex(tm["A"][p - 1], tm["B"][p - 1])) * p ** (1 + be) / dc["C"] for p in (50, 500)]
        assert abs(r[1] - 1) < 0.02 and abs(r[1] - 1) < abs(r[0] - 1)
    assert abs(w.trig_decay(2, 0.5, 0.4)["C"]) < 1e-15


# ------------------------------------------------------------------ line distribution
def test_line_density_and_distribution():
    x = [-4 + 0.1 * i for i in range(81)]
    f = w.dsgn(x, 0, 0.5, 2, 0.5, 0.5)
    assert max(abs(a - math.exp(-v * v) / math.sqrt(math.pi)) for a, v in zip(f, x)) < 1e-15
    tot = simpson(lambda xs: w.dsgn(xs, 0, 1, 1.3, 0.8, 0.3), -60, 60, 200000)
    assert abs(tot - 1) < 1e-9
    assert abs(w.psgn([0.0], 0, 1, 1.3, 0.3, 0.8)[0] - 0.3 / 1.1) < 1e-14


def test_mixture_representation_by_simulation():
    rng = random.Random(101)
    for be, s, l in ((0.3, 0.4, 0.3), (2.5, 1.0, 0.7)):
        x = w.rsgn(20000, 0, s, be, l, 1 - l, rng=rng)
        u = w.psgn(x, 0, s, be, l, 1 - l)
        assert ks_stat(u) < 1.95 / math.sqrt(len(u))
        frac = sum(v < 0 for v in x) / len(x)
        assert abs(frac - l) < 4 * math.sqrt(l * (1 - l) / len(x))


# ------------------------------------------------------------------ wrapped distribution
def test_wrapped_density_matches_image_sums():
    worst = 0.0
    for r in rows("ref_dw.csv"):
        s = 0.4 if r["case"] == "1" else math.sqrt(0.21)
        v = w.dwsgn([float(r["theta"])], 0.7, float(r["beta"]), 0.3, sigma=s)[0]
        worst = max(worst, abs(v / float(r["dens"]) - 1))
    assert worst < 1e-12
    th = [TWO_PI * j / 12 for j in range(13)]
    for be in (0.1, 0.3, 1):
        a = w.dwsgn(th, 0.7, be, 0.3, K=3, sigma=0.4)
        b = w.dwsgn(th, 0.7, be, 0.3, K=400, sigma=0.4)
        assert max(abs(x / y - 1) for x, y in zip(a, b)) < 1e-12
    lg = w.dwsgn(th, 1, 0.5, 0.3, log=True)
    d = w.dwsgn(th, 1, 0.5, 0.3, sigma=math.sqrt(0.21))
    assert max(abs(a - math.log(b)) for a, b in zip(lg, d)) < 1e-13


def test_wrapped_density_integrates_to_one_and_wrapped_normal():
    for be in (0.2, 1, 4):
        tot = simpson(lambda xs: w.dwsgn(xs, 1, be, 0.3, sigma=0.4), 0, TWO_PI, 200000)
        assert abs(tot - 1) < 1e-6
    for t in (0.5, 2.0, 4.5):
        wn = sum(math.exp(-0.5 * ((t + TWO_PI * k) / (math.sqrt(2) * 0.6)) ** 2)
                 for k in range(-80, 81)) / (math.sqrt(2) * 0.6 * math.sqrt(TWO_PI))
        assert abs(w.dwsgn([t], 0, 2, 0.5, sigma=0.6)[0] - wn) < 1e-15


def test_wrapped_distribution_and_quantile():
    for mu, be, lam, s in ((1, 1.5, 0.3, 0.5), (5, 0.2, 0.2, 0.4)):
        q = [0.3, 1.0, 2.5, 4.0, 6.0]
        num = [simpson(lambda xs: w.dwsgn(xs, mu, be, lam, sigma=s), 0, u, 100000) for u in q]
        assert max(abs(a - b) for a, b in zip(w.pwsgn(q, mu, be, lam, sigma=s), num)) < 1e-8
        assert w.pwsgn([-1, 0, TWO_PI, 7], mu, be, lam, sigma=s) == [0, 0, 1, 1]
        p = [0.001, 0.05, 0.3, 0.5, 0.9, 0.999]
        qq = w.qwsgn(p, mu, be, lam, sigma=s)
        assert max(abs(a - b) for a, b in zip(w.pwsgn(qq, mu, be, lam, sigma=s), p)) < 1e-12
    assert abs(w.pwsgn([TWO_PI - 1e-9], 0, 0.3, 0.5)[0] - 1) < 1e-8


def test_wrapped_draws_follow_the_distribution():
    rng = random.Random(23)
    th = w.rwsgn(20000, 4, 0.4, 0.6, sigma=0.3, rng=rng)
    assert all(0 <= v < TWO_PI for v in th)
    assert ks_stat(w.pwsgn(th, 4, 0.4, 0.6, sigma=0.3)) < 1.95 / math.sqrt(len(th))


def test_invalid_parameters():
    with pytest.raises(ValueError):
        w.dwsgn([1.0], 0, -1, 0.5)
    with pytest.raises(ValueError):
        w.dwsgn([1.0], 0, 1, 1.5)
    with pytest.raises(ValueError):
        w.dwsgn([1.0], 0, 1, 0.5, sigma=0)


# ------------------------------------------------------------------ four-parameter likelihood
def test_fit_with_shape_above_one_and_information():
    rng = random.Random(5)
    truth = {"mu": 1.0, "sigma": 0.5, "lam": 0.3, "beta": 1.8}
    th = w.rwsgn(400, 1.0, 1.8, 0.3, sigma=0.5, rng=rng)
    f = w.fit_wsgn(th)
    assert f["se_method"] == "information"
    assert all(abs(f["par"][k] - truth[k]) < 4 * f["se"][k] for k in truth)
    ll = sum(w.dwsgn(th, f["par"]["mu"], f["par"]["beta"], f["par"]["lam"],
                     sigma=f["par"]["sigma"], log=True))
    assert abs(f["logLik"] - ll) < 1e-9
    assert all(abs(math.sqrt(f["vcov"][i][i]) - f["se"][k]) < 1e-12
               for i, k in enumerate(("mu", "sigma", "lam", "beta")))


def vdc(n):
    res = []
    for i in range(1, n + 1):
        s, f, k = 0.0, 0.5, i
        while k > 0:
            s += f * (k % 2)
            k //= 2
            f /= 2
        res.append(s)
    return res


def test_location_is_an_observation_below_shape_one_and_bootstrap():
    rng = random.Random(8)
    th = w.qwsgn(vdc(90), 2, 0.6, 0.4, sigma=0.25)
    with pytest.warns(UserWarning, match="at most one"):
        w.fit_wsgn(th)
    f = w.fit_wsgn(th, se="bootstrap", B=4, rng=rng)
    assert f["par"]["beta"] <= 1
    assert f["par"]["mu"] in [v % TWO_PI for v in th]
    assert f["se_method"] == "bootstrap" and len(f["boot"]) == 4
    assert all(math.isfinite(v) and v >= 0 for v in f["se"].values()) and all(f["se"][k] > 0 for k in ("mu", "sigma", "lam"))


def test_floor_symmetric_profile():
    rng = random.Random(10)
    th = w.rwsgn(150, 1, 0.5, 0.5, sigma=0.05, rng=rng)
    f = w.fit_wsgn(th, beta_floor=1.5, se="none")
    assert f["at_floor"] and abs(f["par"]["beta"] - 1.5) < 1e-5
    fs = w.fit_wsgn(w.rwsgn(200, 1, 2, 0.5, sigma=0.5, rng=rng), symmetric=True)
    assert fs["par"]["lam"] == 0.5 and fs["se"]["lam"] != fs["se"]["lam"]
    th2 = w.rwsgn(60, 1, 1.2, 0.35, sigma=0.4, rng=rng)
    f2 = w.fit_wsgn(th2, se="none")
    pr = w.profile_wsgn(th2, [0.02, f2["par"]["beta"], 3], window=0.4)
    assert abs(pr[1]["logLik"] - f2["logLik"]) < 1e-6
    assert pr[0]["logLik"] > f2["logLik"] > pr[2]["logLik"]


def test_grouped_likelihood_and_fit():
    rng = random.Random(12)
    h = math.pi / 18
    th = [(round(v / h) * h) % TWO_PI for v in w.rwsgn(120, 3, 1.2, 0.4, sigma=0.3, rng=rng)]
    ll = w.grouped_loglik(th, h, 3, 0.3, 0.4, 1.2)
    cnt = [0] * 36
    for v in th:
        cnt[int(round(v / h)) % 36] += 1
    tot = 0.0
    for j in range(36):
        if cnt[j]:
            lo, hi = (j - 0.5) * h, (j + 0.5) * h
            if lo < 0:
                p = w.pwsgn([hi], 3, 1.2, 0.4, sigma=0.3)[0] + 1 - \
                    w.pwsgn([TWO_PI + lo], 3, 1.2, 0.4, sigma=0.3)[0]
            else:
                p = w.pwsgn([hi], 3, 1.2, 0.4, sigma=0.3)[0] - w.pwsgn([lo], 3, 1.2, 0.4, sigma=0.3)[0]
            tot += cnt[j] * math.log(p)
    assert abs(ll - tot) < 1e-9
    g = w.fit_wsgn(th, resolution=h, se="none")
    assert g["grouped"] and ll <= g["logLik"] <= 0
    with pytest.raises(ValueError):
        w.fit_wsgn(th, resolution=0.7)


PERI = [67, 66, 74, 61, 58, 60, 100, 89, 171, 166, 98, 60, 197, 98, 86, 123, 165, 133, 101, 105, 71,
        84, 75, 98, 83, 71, 74, 91, 38, 200, 56]


def test_periwinkle_fit_agrees_with_independent_computation():
    th = [d * math.pi / 180 for d in PERI]
    f = w.fit_wsgn(th, beta_floor=0.3, se="none")
    assert f["par"]["mu"] == 60 * math.pi / 180
    ref = (1.0471975511965979, 0.0810581284965923, 0.12834861700500036, 0.8730181786484532)
    for k, v in zip(("mu", "sigma", "lam", "beta"), ref):
        assert abs(f["par"][k] - v) < 1e-5 * max(1, abs(v))
    assert abs(f["logLik"] + 27.296644722522444) < 1e-7


# ------------------------------------------------------------------ review items
def test_fit_is_the_maximum_over_the_floor_branch():
    th = [d * math.pi / 180 for d in PERI]
    for fl, ll in ((0.10, -21.007222), (0.15, -27.010481)):
        f = w.fit_wsgn(th, beta_floor=fl, se="none")
        assert f["at_floor"]
        assert abs(f["logLik"] - ll) < 1e-6
        assert abs(f["logLik"] - w.profile_wsgn(th, [fl])[0]["logLik"]) < 1e-8
    assert w.fit_wsgn(th, beta_floor=0.10, se="none")["par"]["mu"] == 98 * math.pi / 180
    with pytest.raises(ValueError):
        w.fit_wsgn(th, beta_floor=0.10, se="information")


def test_degenerate_and_non_finite_data_are_refused():
    with pytest.raises(ValueError):
        w.fit_wsgn([2.0] * 10)
    with pytest.raises(ValueError):
        w.fit_wsgn([1.0, 1.0, 2.0, 2.0])
    with pytest.raises(ValueError):
        w.fit_wsgn([0.1, 0.5, math.inf, 2.0], resolution=math.pi / 18, se="none")
    with pytest.raises(ValueError):
        w.grouped_loglik([0.1, 0.5, 1.0, math.inf, 2.0], math.pi / 18, mu=1, sigma=0.4, lam=0.3,
                         beta=0.5)
    with pytest.raises(ValueError):
        w.fit_wsgn([1.0, 1.01] * 5, resolution=math.pi / 18)
    with pytest.raises(ValueError):
        w.dsgn([1.0], beta=-1)


def test_small_shapes_and_lower_tail():
    for be in (0.005, 0.001):
        d = w.dwsgn([0.5, 3.0], 1, be, 0.3, sigma=0.4)
        assert max(abs(v * 2 * math.pi - 1) for v in d) < 1e-12
    g0 = w.dwsgn([0.0], 1, 1.5, 0.3, sigma=0.4)[0]
    for q in (1e-8, 1e-12, 1e-20, 1e-200):
        p = w.pwsgn([q], 1, 1.5, 0.3, sigma=0.4)[0]
        assert abs(p / (q * g0) - 1) < 2 * q + 1e-13
        assert abs(w.qwsgn([p], 1, 1.5, 0.3, sigma=0.4)[0] / q - 1) < 1e-13
    assert abs(w.qwsgn([1e-300], 1, 1.5, 0.3, sigma=0.4)[0] / (1e-300 / g0) - 1) < 1e-12


def test_series_small_shape_not_truncated():
    t = [1.0, 2.0, 3.0]
    zs = w.cfsgn(t, 0, 0.4, 0.05, 0.3, 0.7, method="series")
    zq = w.cfsgn(t, 0, 0.4, 0.05, 0.3, 0.7)
    assert max(abs(a - b) / abs(b) for a, b in zip(zs, zq)) < 1e-11


def test_review_two_items():
    th = [d * math.pi / 180 for d in PERI]
    # a fit on a bound returns no standard errors by default, with a warning saying why
    with pytest.warns(UserWarning, match="bound"):
        f = w.fit_wsgn(th, beta_floor=0.10)
    assert f["se"] is None and f["se_method"] == "none"
    assert all(math.isfinite(v) for v in f["log_rates"].values())
    # the bootstrap draws from the rates, so it runs on such a fit when asked for
    fb = w.fit_wsgn(th, beta_floor=0.10, se="bootstrap", B=2, rng=random.Random(3))
    assert len(fb["boot"]) == 2
    # the profile at the floor equals the fit's floor search
    assert abs(w.profile_wsgn(th, [0.10])[0]["logLik"] - f["logLik"]) < 1e-8
    with pytest.raises(ValueError):
        w.grouped_loglik([1.0, 2.0, 3.0], None, 1, 0.4, 0.3, 1.5)
    for bad in (0, -1, math.inf):
        with pytest.raises(ValueError):
            w.cfsgn([1.0], tol=bad)
    with pytest.raises(ValueError):
        w.trig_moments(2, tol=0)
    with pytest.raises(ValueError):
        w.trig_moments(float("nan"))
    with pytest.warns(UserWarning, match="1e-14"):
        z = w.cfsgn([0.02], 0, 0.5, 50, 0.5, 0.5, tol=1e-16)[0]
    assert abs(z.real - 0.99993468977084714) < 5e-14
    assert w.trig_decay(2, 0.5, 0.4)["C"] == 0.0
    # next to the location at large shapes: p = q g(0) to first order, quantiles relative
    for be in (2, 5, 20, 50):
        g0 = w.dwsgn([0.0], 0, be, 0.5, sigma=1)[0]
        for q in (1e-200, 1e-100, 1e-20, 1e-12):
            p = w.pwsgn([q], 0, be, 0.5, sigma=1)[0]
            assert abs(p / (q * g0) - 1) < 1e-12
            assert abs(w.qwsgn([p], 0, be, 0.5, sigma=1)[0] / q - 1) < 1e-12
    assert w.qwsgn([1.6e-301], 2, 0.3, 0.3, sigma=0.4)[0] > 1e-301
    assert w.psgn([1 + 7e-9], 1, 0.01, 50, 0.5, 0.5)[0] - 0.5 > 1e-7


ONE_SIDED = [
    1.009604364655474, 1.004657904440722, 1.0093145487773232, 1.0256244520063216,
    1.0074606391867311, 0.98580186894609234, 1.0065984637170595, 1.0101988143428431,
    0.99279226427957445, 1.0111255459702329, 0.98472029355692547, 1.0232856837957289,
    1.0288175613974919, 1.0084132950805336, 0.99569417718824671, 0.99498180582504447,
    1.0127538699540894, 0.99464803548308822, 1.0128646185847676, 1.0044964429067433,
    1.0146160194011407, 1.0221974246189103, 1.0216344039323775, 1.0158897549886241,
    1.0181599444614478, 1.0051570129383278, 1.0240024787417437, 0.99100013288575062,
    0.98685336517637889, 1.0325606488507704, 1.0142593392304891, 0.99183099782637141,
    1.0243000988627327, 0.98994829369885839, 1.010084561614935, 0.99538196260318124,
    0.99699785639335869, 1.0030986326568621, 1.0095180398733379, 1.0006461828179507,
    0.99639288371615653, 0.98464278562010277, 0.99642663270803111, 1.0247538535202914,
    1.0082824443277212, 0.99366069768402587, 1.018792802400478, 1.0256692170374269,
    0.98469742137610072, 0.99892417261706734, 0.98298027186974901, 1.0202661035409131,
    0.98884740439255669, 1.013221821395033, 1.0137278000213756, 0.98850451427522712,
    0.98813580326156281, 1.016979420825366, 1.0125863893175369, 0.98804091905557889,
    1.0061193765409069, 1.0023234619731796, 0.99866158143527906, 0.99565836804557883,
    0.99553662614145344, 0.98949500385416889, 0.98890249252086682, 0.9940505156108943,
    1.0137286839473023, 1.0246097090043271, 1.0308027579794943, 0.9942085787616175,
    0.99740814958976787, 0.98273848142108255, 0.98667474134342614, 0.99704928351218436,
    0.99354290159466419, 0.99853528948884573, 1.0144518649693712, 1.0128389619817024,
    1.0146153505757014, 1.000059615498019, 1.0392267092006433, 0.99811343038063516,
    0.98634674769372632, 1.0181248988740339, 1.0046942938763623, 1.0290145269972206,
    1.0009024937052871, 0.98999221129323212, 0.98804536301581991, 1.0361378596259854,
    0.99331013057863404, 1.0005712478484066, 0.98638394106981075, 1.0173503961152328,
    1.0109577183435245, 1.0051770866403977, 1.0328013899255495, 0.99491164143447886
]


def test_one_sided_profile_floor_and_cell_limit():
    # one-sided fit short of the rate bound: flagged, no information standard errors
    # (the data are R's set.seed(3); rwsgn(100, 1, 3, 0.4, sigma = 0.01))
    th = ONE_SIDED
    with pytest.warns(UserWarning, match="one-sided"):
        f = w.fit_wsgn(th)
    assert f["at_one_sided"] and not f["at_rate_bound"]
    assert f["se"] is None and f["se_method"] == "none"
    assert abs(f["logLik"] - 295.6209) < 1e-3
    with pytest.raises(ValueError, match="bound"):
        w.fit_wsgn(th, se="information")
    # the information-matrix check rejects non-finite and indefinite matrices
    nan4 = [float("nan")] * 16
    assert not w._core._info_ok(nan4, nan4, [float("nan")] * 4, 0)
    eye = [1.0 if i % 5 == 0 else 0.0 for i in range(16)]
    assert w._core._info_ok(eye, eye, [1.0] * 4, 0)
    bad = list(eye); bad[5] = -1.0
    assert not w._core._info_ok(bad, bad, [1.0] * 4, 0)
    # profile_wsgn takes its rate bound from beta_floor
    t7 = [1 + k * 1e-14 for k in range(-3, 4)]
    p01 = w.profile_wsgn(t7, [2.0], beta_floor=0.1)[0]["logLik"]
    p1 = w.profile_wsgn(t7, [2.0], beta_floor=1.0)[0]["logLik"]
    assert p01 > 211 and p01 > p1 + 5
    # more than 1e6 cells is refused
    for res in (2 * math.pi / 3e9, 2 * math.pi / 1.5e6):
        with pytest.raises(ValueError, match="1e6"):
            w.grouped_loglik([1.0, 2.0, 3.0], res, 1, 0.4, 0.3, 1.5)
        with pytest.raises(ValueError, match="1e6"):
            w.fit_wsgn([1.0, 2.0, 3.0, 4.0], resolution=res)
    assert math.isfinite(w.grouped_loglik([1.0, 2.0, 3.0], 2 * math.pi / 1e6, 1, 0.4, 0.3, 1.5))


# R's set.seed(204); th <- rwsgn(60, 1, 3, 0.5, sigma = 0.3); th <- 1 + abs(th - 1), in 5-degree cells
GROUPED_ONE_SIDED = [12,12,12,13,14,12,14,18,19,12,15,14,12,16,15,12,18,15,13,18,19,12,13,13,16,13,16,18,12,15,12,16,14,18,16,17,12,17,16,12,16,13,16,19,12,12,18,13,16,15,17,20,16,12,15,12,14,14,12,14]
# R's set.seed(101); th <- rwsgn(40, 1, 3, 0.5, sigma = 0.01); th <- 1 + abs(th - 1)
NEAR_UPPER = [1.0054680894204338,1.0048774145689501,1.0030786355949655,1.0150507290598096,1.0220601381232746,1.0022964678208421,1.0086348537976861,1.0081462251475468,1.0019952876912959,1.0006559014109291,1.007903952935352,1.0137333176826462,1.0052598175413889,1.0029603418235382,1.0104514319588911,1.005400253619396,1.0027982828461215,1.0002755644451882,1.0010531619690808,1.0208086240422118,1.0011417002688394,1.0072434716619632,1.0212476927399723,1.0150696817723264,1.0201813806743609,1.0056492946280895,1.0100931679107685,1.0167770766156776,1.0009008828014525,1.0066337636363361,1.0204885884942336,1.0030627439775626,1.0230527832369389,1.0022468759923446,1.0028363916424894,1.0153042664567189,1.0200423998543617,1.0039017498288949,1.0037264319408346,1.0175036014638863]


def test_observation_within_rounding_of_the_location():
    # 1 - 2^-53 and 1 + 2^-52 are the doubles next to 1; with rates near 5e19 the density a
    # rounding step from the location is zero on either side, never the value at the location
    for x in (1 - 2.0 ** -53, 1 + 2.0 ** -52):
        assert w.dwsgn([x], 1.0, 2.0, 0.5, sigma=1e-20, log=True)[0] < -1e6
    assert w.dwsgn([1.0], 1.0, 2.0, 0.5, sigma=1e-20, log=True)[0] > 40
    t7 = [1 + k * 1e-14 for k in range(-3, 4)]
    assert abs(w.profile_wsgn(t7, [2.0], beta_floor=0.1)[0]["logLik"] - 211.601) < 1e-3
    f = w.fit_wsgn(t7, beta_floor=0.1, se="none")
    pf = w.profile_wsgn(t7, [f["par"]["beta"]], beta_floor=0.1)[0]["logLik"]
    assert abs(pf - f["logLik"]) < 1e-6
    assert abs(f["logLik"] - 214.0974) < 1e-3


def test_grouped_one_sided_fit_is_flagged():
    h = 5 * math.pi / 180
    th = [(k * h) % TWO_PI for k in GROUPED_ONE_SIDED]
    with pytest.warns(UserWarning, match="one-sided") as rec:
        f = w.fit_wsgn(th, resolution=h)
    assert not any("bootstrap" in str(r.message) for r in rec)
    assert f["at_one_sided"] and f["se"] is None
    assert abs(f["log_rates"]["log_a"] - f["log_rates"]["log_b"]) > math.log(1e12)
    with pytest.raises(ValueError, match="one-sided"):
        w.fit_wsgn(th, resolution=h, se="information")


def test_ngrid_cap_symmetric_nan_and_at_upper():
    with pytest.raises(ValueError, match="1e8"):
        w.fit_wsgn([1.0, 2.0, 3.0, 4.0], ngrid=10 ** 9)
    with pytest.raises(ValueError, match="1e8"):
        w.fit_wsgn([1.0, 2.0, 3.0, 4.0], ngrid=2 ** 32 + 200)
    with pytest.warns(UserWarning):
        g = w.fit_wsgn(NEAR_UPPER)
    assert g["par"]["beta"] > 19.99 and g["at_upper"]


# ------------------------------------------------------------------ threads and the cusp polish
def _same_fit(a, b):
    assert a["par"] == b["par"]
    assert a["logLik"] == b["logLik"]
    assert a["log_rates"] == b["log_rates"]
    assert a["nfev"] == b["nfev"]
    for k in ("at_floor", "at_upper", "at_rate_bound", "at_one_sided"):
        assert a[k] == b[k]


def test_fit_and_profile_are_identical_on_one_two_and_four_threads():
    rng = random.Random(32)
    th = w.rwsgn(12, 1.0, 1.6, 0.3, sigma=0.4, rng=rng)
    f1 = w.fit_wsgn(th, se="none", ncores=1)
    _same_fit(f1, w.fit_wsgn(th, se="none", ncores=2))
    _same_fit(f1, w.fit_wsgn(th, se="none", ncores=4))
    p1 = w.profile_wsgn(th, [0.8, 1.6], ncores=1)
    assert p1 == w.profile_wsgn(th, [0.8, 1.6], ncores=2)
    assert p1 == w.profile_wsgn(th, [0.8, 1.6], ncores=4)
    # ties and a grouped fit
    tt = [round(v, 2) for v in w.rwsgn(40, 1.5, 1.2, 0.5, sigma=0.3, rng=rng)]
    _same_fit(w.fit_wsgn(tt, se="none"), w.fit_wsgn(tt, se="none", ncores=2))
    h = 2 * math.pi / 36
    gg = [(round(v / h) * h) % TWO_PI for v in w.rwsgn(100, 3.0, 1.3, 0.45, sigma=0.5, rng=rng)]
    _same_fit(w.fit_wsgn(gg, resolution=h, se="none"),
              w.fit_wsgn(gg, resolution=h, se="none", ncores=2))


def test_chunked_likelihood_evaluations_are_identical_on_one_and_two_threads():
    # n and the number of cells at least 128: the terms are computed in chunks, summed in order
    rng = random.Random(34)
    th = w.rwsgn(130, 1.0, 1.6, 0.3, sigma=0.4, rng=rng)
    p1 = w.profile_wsgn(th, [1.6], ncores=1)
    assert p1 == w.profile_wsgn(th, [1.6], ncores=2)
    assert p1 == w.profile_wsgn(th, [1.6], ncores=4)
    h = 2 * math.pi / 180
    gg = [(round(v / h) * h) % TWO_PI for v in w.rwsgn(60, 3.0, 1.3, 0.45, sigma=0.5, rng=rng)]
    g1 = w.fit_wsgn(gg, resolution=h, se="none")
    _same_fit(g1, w.fit_wsgn(gg, resolution=h, se="none", ncores=2))
    pg = w.profile_wsgn(gg, [1.3], resolution=h, ncores=1)
    assert pg == w.profile_wsgn(gg, [1.3], resolution=h, ncores=2)


def test_ncores_above_256_is_reduced_with_a_warning():
    th = [0.1, 0.5, 0.9, 1.2, 0.3, 0.7]
    with pytest.warns(UserWarning, match="256"):
        p = w.profile_wsgn(th, [1.5], ncores=100000)
    assert p == w.profile_wsgn(th, [1.5])


def test_ncores_is_validated_and_openmp_is_reported():
    assert isinstance(w.HAS_OPENMP, bool)
    th = [0.1, 0.5, 0.9, 1.2, 0.3, 0.7]
    for bad in (0, 1.5, -1, float("nan")):
        with pytest.raises(ValueError, match="ncores"):
            w.fit_wsgn(th, ncores=bad)
        with pytest.raises(ValueError, match="ncores"):
            w.profile_wsgn(th, [1.0], ncores=bad)


def test_one_sided_fit_at_a_cusp_is_polished_to_the_maximum():
    # the seed-101 data of R's tests rounded to 15 significant digits: the final free polish used to
    # collapse along mu and leave the fit 0.0089 below the profile at the upper bound of the shape
    th = [float("%.15g" % v) for v in NEAR_UPPER]
    f = w.fit_wsgn(th, se="none")
    assert abs(f["logLik"] - 150.4029337676) < 1e-8
    assert f["at_upper"] and f["at_one_sided"]
    assert abs(f["par"]["beta"] - 20) < 1e-6
    assert f["logLik"] >= w.profile_wsgn(th, [20.0])[0]["logLik"] - 1e-8
