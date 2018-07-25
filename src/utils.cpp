#include <Rcpp.h>
using namespace Rcpp;

// Global variables -------------------------
double TOL = std::pow(10, -12);

//' Log-sum-exponential trick.
//'
//' @param x A vector to log-sum-exp.
//'
//' @return The log of the sum of the exponential
//'     of the elements in \code{x}.
//'
//' @author David Gerard
// [[Rcpp::export]]
double log_sum_exp(NumericVector x) {
  double max_x = Rcpp::max(x);
  double lse; // the log-sum-exp
  if (max_x == R_NegInf) { // if all -Inf, need to treat this special to avoid -Inf + Inf.
    lse = R_NegInf;
  }
  else {
    lse = max_x + std::log(Rcpp::sum(Rcpp::exp(x - max_x)));
  }
  return lse;
}

//' Log-sum-exponential trick using just two doubles.
//'
//' @param x A double.
//' @param y Another double.
//'
//' @return The log of the sum of the exponential of x and y.
//'
//' @author David Gerard
//'
// [[Rcpp::export]]
double log_sum_exp_2(double x, double y) {
  double z = std::max(x, y);
  double finalval;
  if (z == R_NegInf) {
    finalval = R_NegInf;
  } else {
    finalval = std::log(std::exp(x - z) + std::exp(y - z)) + z;
  }
  return finalval;
}

//' Log of cumulative sum of exponential.
//'
//' @param x A vector to log-cum-sum-exp.
//'
//' @return A NumericVector. Element i is the log of the sum of
//'     the exponential of the first i elements in \code{x}.
//'
//' @author David Gerard
// [[Rcpp::export]]
NumericVector log_cum_sum_exp(NumericVector x) {
  int n = x.length();
  NumericVector y(n);

  y(0) = x(0);
  if (n > 1) {
    for (int i = 1; i < n; i++) {
      y(i) = log_sum_exp_2(y(i - 1), x(i));
    }
  }

  return y;
}

//' The logit function.
//'
//' @param x A double between 0 and 1.
//'
//' @return The logit of \code{x}.
//'
//' @author David Gerard
// [[Rcpp::export]]
double logit(double x) {
  if ((x < TOL) | ((1.0 - x) < TOL)) {
    Rcpp::stop("logit: x must be between 0 and 1.");
  }
  double lv = std::log(x / (1.0 - x));
  return lv;
}

//' The expit (logistic) function.
//'
//' @param x A double.
//'
//' @return The expit (logistic) of \code{x}.
//'
//' @author David Gerard
// [[Rcpp::export]]
double expit(double x) {
  double ev = 1.0 / (1.0 + std::exp(-x));
  return ev;
}
