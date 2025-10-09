// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_RAD_RESULTANT_H
#define WH_RAD_RESULTANT_H

//' Circular resultant vector
//'
//' Computes the circular resultant vector in polar coordinates
//' from a vector of radians. If the resultant length is
//' smaller than \code{DBL_EPSILON}, the mean direction is set to \code{NaN}.
//' @param theta Numeric vector of radians.
//' @return A numeric vector of 2 elements:
//' \itemize{
//'   \item Mean direction in radians
//'   \item Resultant length
//' }
//' @export
// [[Rcpp::export]]
Eigen::VectorXd WH_rad_resultant(const Eigen::VectorXd &theta) {
  double S_sin = theta.array().sin().sum();
  double S_cos = theta.array().cos().sum();
  double mu;
  double R  = std::sqrt(S_sin * S_sin + S_cos * S_cos);
  if (R > DBL_EPSILON) {
    mu = std::atan2(S_sin, S_cos);
  } else {
    mu = NAN;
  }
  Eigen::VectorXd output(2);
  output(0) = mu;
  output(1) = R;
  return output;
}

#endif