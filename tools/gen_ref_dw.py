"""Reference wrapped densities for the wsgnstat tests and vignette (wsgnstat/inst/extdata/ref_dw.csv): each side's image sum by mpmath.nsum
(Euler-Maclaurin extrapolation) at 40 digits; mu = 0.7, lam = 0.3; sigma 0.4 (case 1) or the
reciprocal-rate value sqrt(0.21) (case 2)."""
import sys, math
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import mpmath as mp
from ref_dw import dw
out = open(sys.argv[1], "w")
out.write("beta,case,theta,dens\n")
for be in ("0.1", "0.2", "0.3", "0.5", "1", "2"):
    for case, sig in ((1, mp.mpf("0.4")), (2, mp.sqrt(mp.mpf("0.21")))):
        for th in ("0", "0.7", "1.5", "2.5", "4", "5.5"):
            v = dw(mp.mpf(th), mp.mpf("0.7"), sig, mp.mpf("0.3"), mp.mpf(be))
            out.write(f"{be},{case},{th},{mp.nstr(v, 20)}\n"); out.flush()
