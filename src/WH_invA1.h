// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_INVA1_H
#define WH_INVA1_H

double WH_invA1(double x) {
  double result;
  if (x >= 0.0 && x < 0.53) {
    result = 2.0 * x + std::pow(x, 3) + (5.0 * std::pow(x, 5)) / 6.0;
  } else if (x < 0.85) {
    result = -0.4 + 1.39 * x + 0.43 / (1.0 - x);
  } else {
    result = 1.0 / (std::pow(x, 3) - 4.0 * std::pow(x, 2) + 3.0 * x);
  }
  return result;
}

#endif