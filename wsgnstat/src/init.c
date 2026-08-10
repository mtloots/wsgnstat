/* Native-routine registration for wsgnstat (CRAN-compliant .C interface). */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "wsgn.h"

static const R_CMethodDef CEntries[] = {
  {"wsgn_density",  (DL_FUNC) &wsgn_density,  8},
  {"wsgn_wdensity", (DL_FUNC) &wsgn_wdensity, 7},
  {"wsgn_cf",       (DL_FUNC) &wsgn_cf,       9},
  {"wsgn_trigmom",  (DL_FUNC) &wsgn_trigmom,  6},
  {"wsgn_nll",      (DL_FUNC) &wsgn_nll,      7},
  {NULL, NULL, 0}
};

void R_init_wsgnstat(DllInfo *dll) {
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
