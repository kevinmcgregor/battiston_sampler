// Function to call from R to get Stirling numbers based on Buntine's code

#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include "stable.h"

SEXP stirling(SEXP N, SEXP M, SEXP a) {
  SEXP result = PROTECT(allocVector(REALSXP, 1));
  REAL(result)[0] = asReal(N) + asReal(M);

  SEXP Nint = PROTECT(allocVector(INTSXP, 1));
  SEXP Mint = PROTECT(allocVector(INTSXP, 1));

  stable_t tab = S_make(N, M, N, M, a, S_STABLE);


  UNPROTECT(3);
  return result;
}