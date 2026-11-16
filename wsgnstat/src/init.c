/* Native-routine registration for wsgnstat (CRAN-compliant .C interface). */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "wsgn.h"

static const R_CMethodDef CEntries[] = {
  {"wsgn_density",    (DL_FUNC) &wsgn_density,     8},
  {"wsgn_pdist",      (DL_FUNC) &wsgn_pdist,       8},
  {"wsgn_cf",         (DL_FUNC) &wsgn_cf,         11},
  {"wsgn_trigmom",    (DL_FUNC) &wsgn_trigmom,     8},
  {"wsgn_decay",      (DL_FUNC) &wsgn_decay,       4},
  {"wsgn_wdensity",   (DL_FUNC) &wsgn_wdensity,    9},
  {"wsgn_nll",        (DL_FUNC) &wsgn_nll,         8},
  {"wsgn_wcdf",       (DL_FUNC) &wsgn_wcdf,        8},
  {"wsgn_wquant",     (DL_FUNC) &wsgn_wquant,      8},
  {"wsgn_fit",        (DL_FUNC) &wsgn_fit,        12},
  {"wsgn_profile",    (DL_FUNC) &wsgn_profile,    13},
  {"wsgn_grouped_ll", (DL_FUNC) &wsgn_grouped_ll, 10},
  {"wsgn_fit_batch",  (DL_FUNC) &wsgn_fit_batch,  14},
  {"wsgn_info",       (DL_FUNC) &wsgn_info,        8},
  {NULL, NULL, 0}
};

void R_init_wsgnstat(DllInfo *dll) {
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
