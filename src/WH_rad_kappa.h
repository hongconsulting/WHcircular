// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_RAD_KAPPA_H
#define WH_RAD_KAPPA_H

#include "WH_invA1.h"
#include "WH_rad_mean.h"

//' Circular concentration
//'
//' Computes the von Mises concentration parameter κ of a vector
//' of radians.
//' @param theta Numeric vector of radians.
//' @return The maximum likelihood estimate of κ.
//' @export
// [[Rcpp::export]]
double WH_rad_kappa(const Eigen::VectorXd &theta) {
  double mu = WH_rad_mean(theta);
  if (std::isnan(mu)) return NAN;
  double V = (theta.array() - mu).cos().mean();
  if (V > 0.0) {
    return WH_invA1(V);
  } else {
    return 0.0;
  }
}

#endif