"""Reference wrapped densities: c * sum over both sides of exp(-(rate (d + 2 pi k))^beta), each
side summed by mpmath.nsum (Euler-Maclaurin extrapolation) at 40 digits."""
import mpmath as mp
mp.mp.dps = 40
def dw(theta, mu, sig, lam, be):
    sig = mp.mpf(sig); lam = mp.mpf(lam); be = mp.mpf(be)
    a = lam / sig; b = (1 - lam) / sig
    c = lam * (1 - lam) / (sig * mp.gamma(1 + 1 / be))
    t = mp.fmod(mp.mpf(theta) - mu, 2 * mp.pi)
    if t < 0: t += 2 * mp.pi
    R = mp.nsum(lambda k: mp.exp(-(a * (t + 2 * mp.pi * k)) ** be), [0, mp.inf], method='e')
    Lf = mp.nsum(lambda k: mp.exp(-(b * (2 * mp.pi - t + 2 * mp.pi * k)) ** be), [0, mp.inf], method='e')
    return c * (R + Lf)
