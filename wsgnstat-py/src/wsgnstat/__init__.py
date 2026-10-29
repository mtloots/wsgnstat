"""wsgnstat: the skew generalised normal distribution on the line and its wrapped version on the
circle, on a pure-C back-end shared with the R package 'wsgnstat'. Pure standard library."""
from ._core import (dsgn, psgn, dwsgn, pwsgn, qwsgn, cfsgn, trig_moments, nll, rsgn, rwsgn,
                    fit_wsgn_trig, fit_wsgn_mle)
__all__ = ["dsgn", "psgn", "dwsgn", "pwsgn", "qwsgn", "cfsgn", "trig_moments", "nll", "rsgn",
           "rwsgn", "fit_wsgn_trig", "fit_wsgn_mle"]
__version__ = "0.1.0"
