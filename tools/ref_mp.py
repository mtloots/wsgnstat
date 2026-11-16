"""High-precision reference for G(w) (normalised by Gamma(1+1/beta)): ascending series for beta > 1,
descending for beta < 1, mpmath at 60 digits (the cancellation is absorbed by the working precision)."""
import mpmath as mp
mp.mp.dps = 80
def Gn(w, be):
    w = mp.mpf(w); be = mp.mpf(be)
    if w == 0: return mp.mpc(1)
    s = mp.mpc(0); k = 0
    if be > 1:
        while True:
            term = mp.gamma((k+1)/be) * (1j*w)**k / mp.factorial(k) / be
            s += term; k += 1
            if k > 30 and abs(term) < mp.mpf(10)**-60 * abs(s): break
    elif be < 1:
        while True:
            term = (-1)**k * mp.gamma(k*be+1) * mp.expj(mp.pi*(k*be+1)/2) * w**(-(k*be+1)) / mp.factorial(k)
            s += term; k += 1
            if k > 30 and abs(term) < mp.mpf(10)**-60 * abs(s): break
    else:
        s = 1/(1-1j*w)
    return s / mp.gamma(1 + 1/be)
def cf(t, mu, sig, lam, be):
    a = mp.mpf(lam)/sig; b = (1-mp.mpf(lam))/sig
    return mp.expj(t*mu) * ((1-mp.mpf(lam)) * Gn(t/a, be) + mp.mpf(lam) * mp.conj(Gn(t/b, be)))
