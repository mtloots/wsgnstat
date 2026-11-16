"""Reference characteristic-function values for the wsgnstat tests (wsgnstat/inst/extdata/ref_cf.csv): the
convergent series of G(w) = int_0^inf exp(i w u - u^beta) du (descending for beta < 1, ascending for
beta > 1; 1/(1 - i w) at beta = 1) in mpmath at 120 digits. Case 1: sigma 0.4, lam 0.3; case 2: the
reciprocal-rate member sigma = sqrt(0.21), lam 0.3. Usage: python3 gen_ref_cf.py ref_cf.csv"""
import sys, csv
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import mpmath as mp
import ref_mp
ref_mp.mp.mp.dps = 120
out = open(sys.argv[1], "w")
out.write("beta,case,t,re,im\n")
for be in ("0.1", "0.2", "0.3", "0.4", "0.5", "0.55", "0.6", "0.7", "0.8", "1", "1.5", "2", "3"):
    for case, sig in ((1, mp.mpf("0.4")), (2, mp.sqrt(mp.mpf("0.21")))):
        for t in (1, 2, 3):
            z = ref_mp.cf(t, 0, sig, mp.mpf("0.3"), mp.mpf(be))
            out.write(f"{be},{case},{t},{mp.nstr(z.real, 20)},{mp.nstr(z.imag, 20)}\n")
            out.flush()
