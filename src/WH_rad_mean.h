// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_RAD_MEAN_H
#define WH_RAD_MEAN_H

#include "WH_rad_resultant.h"

//' Circular mean
//'
//' Computes the circular mean of a vector of radians.
//' @param theta Numeric vector of radians.
//' @return The circular mean in radians.
//' @examples
//' y <- c("01:00","01:15","01:30","01:45","02:00",
//'        "23:00","23:15","23:30","23:45","00:00")
//' X <- c(rep(0, 5), rep(1, 5))
//' theta <- WH_HHMM_to_rad(y)
//' print(WH_rad_to_HHMM(WH_rad_sd(theta)))
//' @export
// [[Rcpp::export]]
double WH_rad_mean(const Eigen::VectorXd &theta) {
  return WH_rad_resultant(theta)[0];
}

#endif