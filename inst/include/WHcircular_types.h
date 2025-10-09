// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifndef WH_REGFIT_H
#define WH_REGFIT_H

struct WH_regfit {
  Eigen::MatrixXd vcov;
  Eigen::VectorXd coef;
  Eigen::VectorXd SE;
  Eigen::VectorXd p;
  double n;
  double df;
  double loglik;
  WH_regfit()
    : vcov(), coef(), SE(), p(),
      n(NAN), df(NAN), loglik(NAN) {}
};

struct WH_regfit_list {
private:
  std::vector<WH_regfit> regfits;
public:
  WH_regfit_list(size_t n) : regfits(n) {}
  WH_regfit& operator[](size_t index) {
    return regfits.at(index);
  }
};

namespace Rcpp {
template <> inline SEXP wrap(const WH_regfit& fit) {
  return Rcpp::List::create(
    Rcpp::Named("coef") = fit.coef,
    Rcpp::Named("SE") = fit.SE,
    Rcpp::Named("p") = fit.p,
    Rcpp::Named("n") = fit.n,
    Rcpp::Named("df") = fit.df,
    Rcpp::Named("loglik") = fit.loglik,
    Rcpp::Named("vcov") = fit.vcov
  );
}
template <> inline WH_regfit as(SEXP data) {
  Rcpp::List list(data);
  WH_regfit fit;
  fit.coef = Rcpp::as<Eigen::VectorXd>(list["coef"]);
  fit.SE = Rcpp::as<Eigen::VectorXd>(list["SE"]);
  fit.p = Rcpp::as<Eigen::VectorXd>(list["p"]);
  fit.n = Rcpp::as<double>(list["n"]);
  fit.df = Rcpp::as<double>(list["df"]);
  fit.loglik = Rcpp::as<double>(list["loglik"]);
  fit.vcov = Rcpp::as<Eigen::MatrixXd>(list["vcov"]);
  return fit;
}
template <> inline SEXP wrap(const WH_regfit_list& obj) {
  return Rcpp::XPtr<WH_regfit_list>(new WH_regfit_list(obj), true);
}
template <> inline WH_regfit_list as(SEXP sexp) {
  Rcpp::XPtr<WH_regfit_list> ptr(sexp);
  return *ptr;
}
}

// [[Rcpp::export]]
inline WH_regfit WH_regfit_list_get(Rcpp::XPtr<WH_regfit_list> list, int index) {
  return (*list)[index];
}

// [[Rcpp::export]]
inline void WH_regfit_list_set(Rcpp::XPtr<WH_regfit_list> list, int index, const WH_regfit& regfit) {
  (*list)[index] = regfit;
}

#endif