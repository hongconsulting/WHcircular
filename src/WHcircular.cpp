// [[Rcpp::depends(RcppEigen)]]
#include "WHcircular_types.h"

#include <RcppEigen.h>
#include "WH_distr_chisq.h"
#include "WH_constexpr.h"
#include "WH_invA1.h"
#include "WH_rad_kappa.h"
#include "WH_rad_mean.h"
#include "WH_rad_resultant.h"
// #include "WH_regfit_Rcpp.h"
// #include "WH_regfit.h"
#include "WH_select_bool.h"

double WH_A1(double k) {
  return std::cyl_bessel_i(1, k) / std::cyl_bessel_i(0, k);
}

//' Circular ANOVA
//'
//' One-way analysis of variance for circular data in radians using the
//' Cordeiro–Paula–Botter¹ likelihood-ratio test.
//' @param theta Numeric vector of radians.
//' @param group Integer vector of group identifiers.
//' @return \emph{P}-value of the test statistic.
//' @references
//' 1. Cordeiro, G.M., Paula, G.A. and Botter, D.A., 1994. Improved likelihood
//' ratio tests for dispersion models. \emph{International Statistical Review}, pp. 257–274.
//' @examples
//' y <- c("01:00","01:15","01:30","01:45","02:00",
//'        "23:00","23:15","23:30","23:45","00:00")
//' X <- c(rep(0, 5), rep(1, 5))
//' theta <- WH_HHMM_to_rad(y)
//' print(WH_CordeiroPaulaBotter(theta, X))
//' @export
// [[Rcpp::export]]
double WH_CordeiroPaulaBotter(const Eigen::VectorXd &theta, const Eigen::VectorXi &group) {
  Eigen::VectorXd theta_mod = theta.unaryExpr([](double x){return std::fmod(x, 2.0 * WH::_PI);});
  int n = theta_mod.size();
  std::set<int> u_groups(group.data(), group.data() + n);
  int n_groups = u_groups.size();
  std::vector<int> index(u_groups.begin(), u_groups.end());
  Eigen::VectorXd n_i = Eigen::VectorXd::Zero(n_groups);
  Eigen::VectorXd mu_i = Eigen::VectorXd::Zero(n_groups);
  Eigen::VectorXd R_i = Eigen::VectorXd::Zero(n_groups);
  for (int i = 0; i < n_groups; i++) {
    Eigen::Array<bool, Eigen::Dynamic, 1> mask = (group.array() == index[i]);
    Eigen::VectorXd subset = WH_select_bool(theta_mod, mask);
    n_i(i) = subset.size();
    Eigen::VectorXd resultant = WH_rad_resultant(subset);
    mu_i(i) = resultant[0];
    R_i(i) = resultant[1];
  }
  double mu_all = WH_rad_mean(theta_mod);
  double kappa = WH_rad_kappa(theta_mod);
  double inv_sum = (n_i.cwiseInverse()).sum();
  double stat1 = 1.0 - (1.0 / (4.0 * kappa)) * WH_A1(kappa) * (inv_sum - 1.0 / n);
  double stat2 = (R_i.array() * (1.0 - (mu_i.array() - mu_all).cos())).sum();
  return 1.0 - WH_pchisq(2.0 * kappa * stat1 * stat2, n_groups - 1.0);
}  

//' Convert "HH:MM" time-of-day strings to radians
//'
//' Converts times-of-day in "HH:MM" string format (for example, "14:00" for 2:00 PM)
//' to radians.
//' @param input String vector of times in "HH:MM" format.
//' @return A numeric vector of radians.
//' @export
// [[Rcpp::export]]
Eigen::VectorXd WH_HHMM_to_rad(const std::vector<std::string>& input) {
  int n = input.size();
  Eigen::VectorXd output(n);
  for (int i = 0; i < n; ++i) {
    const std::string& HHMM = input[i];
    size_t pos = HHMM.find(':');
    if (pos == std::string::npos) {
      output(i) = std::numeric_limits<double>::quiet_NaN();
      continue;
    }
    try {
      int HH = std::stoi(HHMM.substr(0, pos));
      int MM = std::stoi(HHMM.substr(pos + 1));
      if (HH < 0 || HH > 23 || MM < 0 || MM > 59) {
        output(i) = std::numeric_limits<double>::quiet_NaN();
      } else {
        int mins = HH * 60 + MM;
        output(i) = (mins / 1440.0) * 2.0 * WH::_PI;
      }
    } catch (...) {
      output(i) = std::numeric_limits<double>::quiet_NaN();
    }
  }
  return output;
}

//' Convert integer times-of-day to radians
//'
//' Converts times-of-day in HHMM integer format (for example, 1400 for 2:00 PM)
//' to radians.
//' @param x Integer vector of times in HHMM format.
//' @return A numeric vector of radians.
//' @export
// [[Rcpp::export]]
Eigen::VectorXd WH_int_to_rad(const Eigen::VectorXi &x) {
  int n = x.size();
  Eigen::VectorXd output(n);
  for (int i = 0; i < n; ++i) {
    int HH = x(i) / 100;
    int MM = x(i) % 100;
    int mins = HH * 60 + MM;
    output(i) = (mins / 1440.0) * 2.0 * WH::_PI;
  }
  return output;
}

//' Convert radians to signed "HH:MM" time differences
//'
//' Converts radians to time-of-day differences in "HH:MM" string format,
//' rounded to the nearest minute. Negative values are prefixed with a minus sign.
//' @param input Numeric vector of radians.
//' @return A string vector of signed time-of-day differences in "HH:MM" format.
//' @export
// [[Rcpp::export]]
std::vector<std::string> WH_rad_to_dHHMM(const Eigen::VectorXd &input) {
  int n = input.size();
  std::vector<std::string> output(n);
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(input(i))) {
      output[i] = "NaN";
      continue;
    }
    double M = input(i) / (2.0 * WH::_PI) * 1440.0;
    std::string sign;
    if (M < 0) sign = "\u2212"; else sign = "";
    int absM = std::round(std::abs(M));
    int absHH = absM / 60;
    int absMM = absM % 60;
    std::string sHH = std::to_string(absHH);
    if (absHH < 10) sHH = "0" + sHH;
    std::string sMM = std::to_string(absMM);
    if (absMM < 10) sMM = "0" + sMM;
    output[i] = sign + sHH + ":" + sMM;
  }
  return output;
}
 
//' Convert radians to "HH:MM" times
//'
//' Converts radians to times-of-day in "HH:MM" string format,
//' rounded to the nearest minute.
//' @param input Numeric vector of radians.
//' @return A string vector of times-of-day in "HH:MM" format.
//' @export
// [[Rcpp::export]]
std::vector<std::string> WH_rad_to_HHMM(const Eigen::VectorXd &input) {
  int n = input.size();
  std::vector<std::string> output(n);
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(input(i))) {
      output[i] = "NaN";
      continue;
    }
    double M = std::fmod(input(i), 2.0 * WH::_PI);
    if (M < 0) M += 2.0 * WH::_PI;
    M = (M / (2.0 * WH::_PI)) * 1440.0;
    int roundM = std::round(M);
    if (roundM == 1440) roundM = 0; // wrap midnight
    int HH = roundM / 60;
    int MM = roundM % 60;
    std::string sHH = std::to_string(HH);
    if (HH < 10) sHH = "0" + sHH;
    std::string sMM = std::to_string(MM);
    if (MM < 10) sMM = "0" + sMM;
    output[i] = sHH + ":" + sMM;
  }
  return output;
}

//' Circular-linear regression
//'
//' Fits a circular-linear regression¹ model for outcomes in radians.
//' @param y Numeric vector of radians.
//' @param X Numeric matrix of predictors.
//' @param tol Convergence tolerance. Default = 10⁻¹⁰.
//' @return A \code{WH_regfit} list with elements:
//' \itemize{
//'   \item \code{coef}: vector of estimated coefficients, ordered as
//'         κ, μ, \bold{β}.
//'   \item \code{SE}: vector of estimated standard errors, ordered as
//'         κ, μ, \bold{β}.
//'   \item \code{p}: \emph{p}-values of the betas.
//'   \item \code{n}: sample size.
//'   \item \code{loglik}: maximized log-likelihood.
//'   \item \code{vcov}: variance–covariance matrix of the betas.
//' }
//' @references
//' 1. Fisher, N.I. and Lee, A.J., 1992. Regression models for an angular response.
//' \emph{Biometrics}, pp. 665–677.
//' @examples
//' y <- c("01:00","01:15","01:30","01:45","02:00",
//'        "23:00","23:15","23:30","23:45","00:00")
//' X <- c(rep(0, 5), rep(1, 5))
//' theta <- WH_HHMM_to_rad(y)
//' fit <- WH_reg_circlinear(theta, as.matrix(X))
//' print(fit$coef)
//' print(fit$SE)
//' @export
// [[Rcpp::export]]
WH_regfit WH_reg_circlinear(const Eigen::VectorXd &y, const Eigen::MatrixXd &X, double tol = 1e-10) {
  int n = y.size();
  int n_beta = X.cols();
  Eigen::VectorXd y_mod = y.unaryExpr([](double x){return std::fmod(x, 2.0 * WH::_PI);});
  Eigen::VectorXd beta_prev = Eigen::VectorXd::Zero(n_beta);
  double S = (y_mod.array() - (2 * (X * beta_prev).array().atan())).sin().sum() / n;
  double C = (y_mod.array() - (2 * (X * beta_prev).array().atan())).cos().sum() / n;
  double R = std::sqrt(S*S + C*C);
  double mu = std::atan2(S, C);
  double k = WH_invA1(R);
  double diff = tol + 1.0;
  Eigen::VectorXd beta_new;
  double A;
  Eigen::MatrixXd D(n, n_beta);
  while (diff > tol) {
    Eigen::VectorXd u = (y_mod.array() - mu - 2 * (X * beta_prev).array().atan()).sin() * k;
    A = WH_A1(k) * k;
    Eigen::VectorXd g(n);
    for (int i = 0; i < n; i++) {
      double denom = 1.0 + std::pow(beta_prev.dot(X.row(i)), 2);
      g(i) = 2.0 / denom;
    }
    for (int i = 0; i < n; i++) {
      D.row(i) = g(i) * X.row(i);
    }
    Eigen::MatrixXd LHS = A * (D.transpose() * D);
    Eigen::VectorXd RHS = D.transpose() * (u + A * (D * beta_prev));
    beta_new = LHS.ldlt().solve(RHS);
    diff = (beta_new - beta_prev).cwiseAbs().maxCoeff();
    beta_prev = beta_new;
    S = (y_mod.array() - (2 * (X * beta_prev).array().atan())).sin().sum() / n;
    C = (y_mod.array() - (2 * (X * beta_prev).array().atan())).cos().sum() / n;
    R = std::sqrt(S*S + C*C);
    mu = std::atan2(S, C);
    k = WH_invA1(R);
  }
  double loglik = -n * std::log(std::cyl_bessel_i(0, k)) +
    k * (y_mod.array() - mu - 2*(X * beta_new).array().atan()).cos().sum();
  Eigen::MatrixXd vcov = (A * (D.transpose() * D)).inverse();
  Eigen::VectorXd beta_SE = vcov.diagonal().array().sqrt();
  double k_SE = std::sqrt(1.0 / (n * (1.0 - std::pow(WH_A1(k), 2) - WH_A1(k)/k)));
  double mu_SE = 1.0 / std::sqrt((n - n_beta) * k * WH_A1(k));
  Eigen::VectorXd t = beta_new.cwiseQuotient(beta_SE).cwiseAbs();
  Eigen::VectorXd p(t.size());
  for (int i = 0; i < n_beta; i++) {
    p(i) = 1.0 - 0.5 * std::erfc(-t(i) * WH::_SQRT1_2);
  }
  WH_regfit output;
  output.n = n;
  output.coef.resize(2 + n_beta);
  output.coef(0) = k;
  output.coef(1) = mu;
  output.coef.tail(n_beta) = beta_new;
  output.SE.resize(2 + n_beta);
  output.SE(0) = k_SE;
  output.SE(1) = mu_SE;
  output.SE.tail(n_beta) = beta_SE;
  output.p = p;
  output.vcov = vcov;
  output.loglik = loglik;
  return output;
}

//' Circular standard deviation
//'
//' Computes the circular standard deviation of a vector of radians.
//' @param theta Numeric vector of radians.
//' @return The circular standard deviation in radians.
//' @examples
//' y <- c("01:00","01:15","01:30","01:45","02:00",
//'        "23:00","23:15","23:30","23:45","00:00")
//' X <- c(rep(0, 5), rep(1, 5))
//' theta <- WH_HHMM_to_rad(y)
//' print(WH_rad_to_HHMM(WH_rad_sd(theta)))
//' @export
// [[Rcpp::export]]
double WH_rad_sd(const Eigen::VectorXd &theta) {
  double S_sin = theta.array().sin().sum();
  double S_cos = theta.array().cos().sum();
  double rho = std::sqrt(S_sin * S_sin + S_cos * S_cos) / theta.size();
  return std::sqrt(-2 * std::log(rho));
}

