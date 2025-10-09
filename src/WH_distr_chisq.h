// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_DISTR_CHISQ_H
#define WH_DISTR_CHISQ_H

#include "WH_Cephes_SciPy_igam.h"

double WH_pchisq(double x, double df) {
  if( (x < 0.0) || (df < 1.0) ) {
    throw std::invalid_argument("[WH_pchisq] df < 1 or x < 0");
  }
  return(Cephes_SciPy::igam(df/2.0, x/2.0));
}

#endif