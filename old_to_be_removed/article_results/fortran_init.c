#ifndef R_IMR_H
#define R_IMR_H
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <R_ext/RS.h> // F77_SUB
#define _(String) (String)

#define FDEF(name) {#name, (DL_FUNC) &F77_SUB(name),           \
   (int)(sizeof(name##_t)/sizeof(name##_t[0])), name##_t}

void F77_SUB(pcrossprod)(
    int    *nrow,
    int    *ncol,
    int    *nrank,
    double *u,
    double *v,
    int    *irow,
    int    *pcol,
    int    *nomega,
    double *r

);

void F77_SUB(pcrossprodt)(
    int    *nrow,
    int    *ncol,
    int    *nrank,
    double *u,
    double *v,
    int    *irow,
    int    *pcol,
    int    *nomega,
    double *r

);


static R_NativePrimitiveArgType pcrossprod_t[] = {
  INTSXP,  /* nrow   */
  INTSXP,  /* ncol   */
  INTSXP,  /* nrank  */
  REALSXP, /* u      */
  REALSXP, /* v      */
  INTSXP,  /* irow   */
  INTSXP,  /* pcol   */
  INTSXP,  /* nomega */
  REALSXP  /* r      */
};


static R_NativePrimitiveArgType pcrossprodt_t[] = {
  INTSXP,  /* nrow   */
  INTSXP,  /* ncol   */
  INTSXP,  /* nrank  */
  REALSXP, /* u      */
  REALSXP, /* v      */
  INTSXP,  /* irow   */
  INTSXP,  /* pcol   */
  INTSXP,  /* nomega */
  REALSXP  /* r      */
};

static R_FortranMethodDef fMethods[] = {
  FDEF(pcrossprod),
  FDEF(pcrossprodt),
  {NULL, NULL, 0}
};

#ifdef __GNUC__
extern void Rcpp_R_init_IMR(DllInfo *dll) __attribute__((weak));
#else
extern void Rcpp_R_init_IMR(DllInfo *dll);
#endif

/* Undo the macro (if any) so our function name stays R_init_IMR */
#ifdef R_init_IMR
#undef R_init_IMR
#endif


void R_init_IMR(DllInfo *dll){
  if(Rcpp_R_init_IMR) Rcpp_R_init_IMR(dll);
  R_registerRoutines(dll, NULL, NULL, fMethods, NULL);
  R_useDynamicSymbols(dll, FALSE);
}


/* .Call calls */
// generated with tools::package_native_routine_registration_skeleton(".")



#endif
