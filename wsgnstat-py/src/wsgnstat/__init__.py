"""wsgnstat: the two-piece exponential power distribution on the line and its wrapped version on
the circle, treated as a uniform scale mixture with generalised gamma mixing, on a pure-C back-end
shared with the R package 'wsgnstat'. The name follows the authors' earlier usage, the (wrapped)
skew generalised normal distribution. Pure standard library."""
from ._core import (dsgn, psgn, cfsgn, rsgn, dwsgn, pwsgn, qwsgn, rwsgn, trig_moments,
                    trig_decay, nll, fit_wsgn, profile_wsgn, grouped_loglik, fit_wsgn_trig,
                    fit_wsgn_mle, HAS_OPENMP)
__all__ = ["dsgn", "psgn", "cfsgn", "rsgn", "dwsgn", "pwsgn", "qwsgn", "rwsgn", "trig_moments",
           "trig_decay", "nll", "fit_wsgn", "profile_wsgn", "grouped_loglik", "fit_wsgn_trig",
           "fit_wsgn_mle", "HAS_OPENMP"]
__version__ = "0.1.1"
