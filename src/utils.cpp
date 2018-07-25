#include <RcppArmadillo.h>
#include <updog.h> // for log_sum_exp_2
using namespace Rcpp;

// Global variables -------------------------
double TOL = std::pow(10, -12);

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
      y(i) = updog::log_sum_exp_2(y(i - 1), x(i));
    }
  }

  return y;
}
