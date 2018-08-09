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


//' Calculate sum_n sum_j sum_l cross_mat(j, l) * pm1(n, j) * pm2(n, l)
//'
//' @param pm1 A matrix of posterior probabilities for SNP 1. Element (i, j)
//'     is the posterior probability of genotype j - 1 in individual i.
//' @param pm2 A matrix of posterior probabilities for SNP 2. Element (i, j)
//'     is the posterior probability of genotype j - 1 in individual i.
//' @param crossmat A matrix of probabilities. Element (j,k) is the
//'     prior probability of genotype j-1 in SNP 1 and k-1 in SNP 2.
//'
//' @author David Gerard
// [[Rcpp::export]]
double sum_out_prods(const NumericMatrix& pm1,
                     const NumericMatrix& pm2,
                     const NumericMatrix& crossmat) {
  // check input -------------------------------
  if (pm1.nrow() != pm2.nrow()) {
    Rcpp::stop("sum_out_prods: pm1 and pm2 should have the same dimension");
  }
  if (pm1.ncol() != pm2.ncol()) {
    Rcpp::stop("sum_out_prods: pm1 and pm2 should have the same dimension");
  }
  if (pm1.ncol() != crossmat.ncol()) {
    Rcpp::stop("sum_out_prods: pm1 and crossmat should have the same number of columns");
  }
  if (crossmat.ncol() != crossmat.nrow()) {
    Rcpp::stop("sum_out_prods: crossmat should have the same number of rows as columns.");
  }

  double cov = 0.0;
  int n = pm1.nrow();
  int ploidy = pm1.ncol() - 1;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= ploidy; j++) {
      for (int ell = 0; ell <= ploidy; ell++) {
        cov += crossmat(j, ell) * pm1(i, j) * pm2(i, ell);
      }
    }
  }

  return cov;
}







