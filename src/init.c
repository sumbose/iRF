#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
   */

/* .C calls */
extern void classForest(void *, void *, void *, void *, void *, void *, void *,
    void *, void *, void *, void *, void *, void *, void *, void *, void *, 
    void *, void *, void *, void *, void *, void *, void *, void *);
extern void classRF(void *, void *, void *, void *, void *, void *, void *, 
    void *, void *, void *, void *, void *, void *, void *, void *, void *, 
    void *, void *, void *, void *, void *, void *, void *, void *, void *, 
    void *, void *, void *, void *, void *, void *, void *, void *, void *, 
    void *, void *, void *, void *, void *, void *, void *, void *, void *, 
    void *, void *);
extern void regForest(void *, void *, void *, void *, void *, void *, void *, 
    void *, void *, void *, void *, void *, void *, void *, void *, void *, 
    void *, void *, void *, void *, void *);
extern void regRF(void *, void *, void *, void *, void *, void *, void *, void *, 
    void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, 
    void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, 
    void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, 
    void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP iRF_RIT_1class(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP iRF_RIT_2class(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
  {"classForest", (DL_FUNC) &classForest, 24},
  {"classRF",     (DL_FUNC) &classRF,     45},
  {"regForest",   (DL_FUNC) &regForest,   21},
  {"regRF",       (DL_FUNC) &regRF,       46},
  {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
  {"iRF_RIT_1class", (DL_FUNC) &iRF_RIT_1class,  9},
  {"iRF_RIT_2class", (DL_FUNC) &iRF_RIT_2class, 11},
  {NULL, NULL, 0}
};

void R_init_iRF(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
