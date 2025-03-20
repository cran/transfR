#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern void gdist(double *coord1, double *coord2, int n1, int n2,
                  int proj, int rescale, int diag, int nthreads, double *mdist);
SEXP c_gdist(SEXP coord1, SEXP coord2, SEXP proj, SEXP rescale, SEXP diag, SEXP nthreads) {
    SEXP dim1 = PROTECT(getAttrib(coord1, R_DimSymbol));
    SEXP dim2 = PROTECT(getAttrib(coord2, R_DimSymbol));
    const int n1 = INTEGER(dim1)[0];
    const int n2 = INTEGER(dim2)[0];
    SEXP mdist = PROTECT(allocVector(REALSXP, 1));
    gdist(REAL(coord1), REAL(coord2), n1, n2, asInteger(proj), asInteger(rescale), asInteger(diag), asInteger(nthreads), REAL(mdist));
    UNPROTECT(3);
    return mdist;
}

extern void similarity(double *Rn, int nt, int nb, int crit, int nthreads, double *sim_matrix);
SEXP c_similarity(SEXP Rn, SEXP crit, SEXP nthreads) {
    SEXP dim = PROTECT(getAttrib(Rn, R_DimSymbol));
    const int nt = INTEGER(dim)[0];
    const int nb = INTEGER(dim)[1];
    SEXP sim_matrix = PROTECT(allocMatrix(REALSXP, nb, nb));
    similarity(REAL(Rn), nt, nb, asInteger(crit), asInteger(nthreads), REAL(sim_matrix));
    UNPROTECT(2);
    return sim_matrix;
}

extern void convolution(double *rn, double *uh, int nrn, int nuh, double *q);
SEXP c_convolution(SEXP rn, SEXP uh) {
    const int nrn = length(rn);
    const int nuh = length(uh);
    SEXP q = PROTECT(allocVector(REALSXP, nrn + nuh));   
    convolution(REAL(rn), REAL(uh), nrn, nuh, REAL(q));
    UNPROTECT(1);
    return q;
}

static const R_CallMethodDef CallEntries[] = {
    {"c_gdist", (DL_FUNC) &c_gdist, 6},
    {"c_similarity", (DL_FUNC) &c_similarity, 3},
    {"c_convolution", (DL_FUNC) &c_convolution, 2},
    {NULL, NULL, 0}
};

void R_init_transfR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
