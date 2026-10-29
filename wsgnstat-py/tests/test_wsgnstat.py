"""Tests of the Python front. They defend the same correctness claims as the R tests: the densities
integrate to one, the circular law reduces to the wrapped normal at beta = 2 and lam = 0.5, the
distribution functions are the integrals of the densities, the quantile function inverts the
distribution function, the sampler follows the law, and the estimators recover the parameters that
generated the data."""
import math
import random
import warnings

import pytest

import wsgnstat as w

TWO_PI = 2.0 * math.pi


def simpson(f, a, b, n=20000):
    """Composite Simpson on [a, b] with n (even) panels, f vectorised over a list."""
    h = (b - a) / n
    xs = [a + h * i for i in range(n + 1)]
    ys = f(xs)
    s = ys[0] + ys[-1] + 4 * sum(ys[1:-1:2]) + 2 * sum(ys[2:-1:2])
    return s * h / 3.0


def norm_pdf(x, s):
    return math.exp(-0.5 * (x / s) ** 2) / (s * math.sqrt(2 * math.pi))


def norm_cdf(x, s):
    return 0.5 * math.erfc(-x / (s * math.sqrt(2)))


def wrapped_normal_pdf(t, s, k=60):
    return sum(norm_pdf(t + TWO_PI * j, s) for j in range(-k, k + 1))


@pytest.mark.parametrize("par", [(0, 1, 2, 0.5, 0.5), (0, 1, 1.3, 0.8, 0.3), (1, 2, 3, 0.2, 0.9)])
def test_line_density_integrates_to_one(par):
    tot = simpson(lambda xs: w.dsgn(xs, *par), -40, 40, 80000)
    assert abs(tot - 1) < 1e-4
    assert min(w.dsgn([-5, 0, 5], *par)) >= 0


@pytest.mark.parametrize("par", [(0, 2, 0.5), (1, 1.5, 0.3), (0, 3, 0.7)])
def test_circular_density_integrates_to_one(par):
    tot = simpson(lambda ts: w.dwsgn(ts, *par), 0, TWO_PI, 20000)
    assert abs(tot - 1) < 1e-6


def test_beta2_is_the_wrapped_normal():
    th = [0.5, 2.0, 4.5, 6.0]
    got = w.dwsgn(th, 0, 2, 0.5)
    for g, t in zip(got, th):
        assert abs(g - wrapped_normal_pdf(t, math.sqrt(0.5))) < 1e-9


def test_trig_moments_closed_form():
    tm = w.trig_moments(4, 0, 2, 0.5)
    for p in range(1, 5):
        assert abs(tm["A"][p - 1] - math.exp(-p * p / 4)) < 1e-9
        assert abs(tm["B"][p - 1]) < 1e-9


def test_characteristic_function():
    assert abs(w.cfsgn([0], 0, 1, 2, 0.5, 0.5)[0] - 1) < 1e-10
    # Fourier transform of the density at beta = 1.3, a skew member
    for t in (0.4, 1.1):
        re = simpson(lambda xs: [math.cos(t * x) * d for x, d in
                                 zip(xs, w.dsgn(xs, 0, 1, 1.3, 0.8, 0.3))], -25, 25, 100000)
        im = simpson(lambda xs: [math.sin(t * x) * d for x, d in
                                 zip(xs, w.dsgn(xs, 0, 1, 1.3, 0.8, 0.3))], -25, 25, 100000)
        cf = w.cfsgn([t], 0, 1, 1.3, 0.8, 0.3)[0]
        assert abs(cf - complex(re, im)) < 1e-4


def test_line_distribution_function():
    par = (0, 1, 1.3, 0.3, 0.8)
    for q in (-3, -0.5, 0.4, 2.5):
        num = simpson(lambda xs: w.dsgn(xs, *par), -60, q, 200000)
        assert abs(w.psgn([q], *par)[0] - num) < 1e-7
    assert abs(w.psgn([0], *par)[0] - 0.3 / 1.1) < 1e-14
    for q in (-1.2, 0.4, 2):
        assert abs(w.psgn([q], 0, 0.5, 2, 0.5, 0.5)[0] - norm_cdf(q, math.sqrt(0.5))) < 1e-14
    assert w.psgn([-math.inf, math.inf], *par) == [0.0, 1.0]


def test_wrapped_distribution_function():
    for par in ((1, 1.5, 0.3), (5, 0.8, 0.2), (3, 8, 0.9)):
        for q in (0.3, 2.5, 6.0):
            num = simpson(lambda ts: w.dwsgn(ts, *par), 0, q, 40000)
            assert abs(w.pwsgn([q], *par)[0] - num) < 1e-7
        assert abs(w.pwsgn([TWO_PI], *par)[0] - 1) < 1e-12
    def wn_cdf(t, s=math.sqrt(0.5)):
        return sum(norm_cdf(t + TWO_PI * j, s) - norm_cdf(TWO_PI * j, s) for j in range(-60, 61))
    for q in (1, 3, 5):
        assert abs(w.pwsgn([q], 0, 2, 0.5)[0] - wn_cdf(q)) < 1e-12


def test_quantile_inverts():
    p = [0.001, 0.05, 0.3, 0.5, 0.9, 0.999]
    for par in ((1, 1.5, 0.3), (5, 0.8, 0.2), (0, 2, 0.5)):
        q = w.qwsgn(p, *par)
        back = w.pwsgn(q, *par)
        assert max(abs(a - b) for a, b in zip(back, p)) < 1e-12
    assert w.qwsgn([0, 1], 1, 1.5, 0.3) == [0.0, TWO_PI]


def test_truncation_is_reported():
    tot = w.pwsgn([TWO_PI], 0, 0.3, 0.5)[0]
    assert tot < 0.95
    assert w.pwsgn([TWO_PI], 0, 0.3, 0.5, K=400)[0] > tot
    with pytest.warns(UserWarning, match="raise K"):
        w.qwsgn([0.99], 0, 0.3, 0.5)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        w.qwsgn([0.5], 0, 0.3, 0.5)


def test_draws_follow_the_law():
    rng = random.Random(23)
    n = 20000
    th = w.rwsgn(n, 1, 1.5, 0.3, rng=rng)
    assert all(0 <= t < TWO_PI for t in th)
    u = sorted(w.pwsgn(th, 1, 1.5, 0.3))
    ks = max(max((i + 1) / n - v, v - i / n) for i, v in enumerate(u))
    assert ks < 1.95 / math.sqrt(n)
    x = w.rsgn(n, 0, 1, 1.3, 0.3, 0.8, rng=rng)
    frac = sum(1 for v in x if v < 0) / n
    assert abs(frac - 0.3 / 1.1) < 4 * math.sqrt(0.25 / n)


@pytest.mark.parametrize("truth", [(1, 1.5, 0.3), (4, 0.8, 0.6), (2, 3, 0.5)])
def test_estimators_recover_parameters(truth):
    rng = random.Random(17)
    th = w.rwsgn(2000, *truth, rng=rng)
    ft = w.fit_wsgn_trig(th)
    names = ("mu", "beta", "lam")
    for k, t in zip(names, truth):
        assert ft["se"][k] > 0
        assert abs(ft["par"][k] - t) < 4 * ft["se"][k]
    fm = w.fit_wsgn_mle(th)
    assert abs(fm["par"]["mu"] - truth[0]) < 0.1
    assert abs(fm["par"]["beta"] - truth[1]) / truth[1] < 0.15   # MC sd 3.5 per cent at beta = 3, n = 2000
    assert abs(fm["par"]["lam"] - truth[2]) < 0.06
    p = ft["par"]
    assert fm["logLik"] >= -w.nll(th, p["mu"], p["beta"], p["lam"]) - 1e-8


def test_bootstrap_standard_errors():
    rng = random.Random(31)
    th = w.rwsgn(400, 1, 1.5, 0.3, rng=rng)
    fm = w.fit_wsgn_mle(th, B=15, rng=rng)
    ft = w.fit_wsgn_trig(th)
    for k in ("mu", "beta", "lam"):
        assert 0 < fm["se"][k]
        assert ft["se"][k] / 5 < fm["se"][k] < 5 * ft["se"][k]
