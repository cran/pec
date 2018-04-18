#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void ccr(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void cindexSRC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void pecCR(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void pec_noinf(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void pec_noinfCR(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void pecResiduals(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void pecResidualsCR(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void pecSRC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void SNull(void *, void *, void *, void *, void *, void *);
extern void survest_cox_aalen(void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"ccr",               (DL_FUNC) &ccr,               19},
    {"cindexSRC",         (DL_FUNC) &cindexSRC,         16},
    {"pecCR",             (DL_FUNC) &pecCR,             12},
    {"pec_noinf",         (DL_FUNC) &pec_noinf,         11},
    {"pec_noinfCR",       (DL_FUNC) &pec_noinfCR,       12},
    {"pecResiduals",      (DL_FUNC) &pecResiduals,      12},
    {"pecResidualsCR",    (DL_FUNC) &pecResidualsCR,    13},
    {"pecSRC",            (DL_FUNC) &pecSRC,            11},
    {"SNull",             (DL_FUNC) &SNull,              6},
    {"survest_cox_aalen", (DL_FUNC) &survest_cox_aalen,  6},
    {NULL, NULL, 0}
};

void R_init_pec(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
