// Functions for use in correlation estimation

#include <RcppArmadillo.h>
#include <math.h> // for M_PI
#include <updog.h> // for log_sum_exp_2
using namespace Rcpp;

// External functions ----------------------
NumericVector log_cum_sum_exp(NumericVector x); // From utils.cpp
double pnormcop(double x, double y, double rho); // From gauss_cop.cpp

// Global variables -------------------------
extern double TOL; // defined in utils.cpp

//' Likelihood of the correlation given the count data.
//'
//' @param atanh_rho The inverse hyperbolic tangent of the correlation.
//'     This is also known as the Fisher Z-transformation of the correlation.
//' @param lX A matrix with the genotype log-likelihoods for the first SNP. Element
//'     (\code{i}, \code{j}) is the log-likelihood
//'     of the data for individual \code{i} given dosage \code{j} for SNP 1.
//' @param lY A matrix with the genotype log-likelihoods for the second SNP. Element
//'     (\code{i}, \code{j}) is the log-likelihood
//'     of the data for individual \code{i} given dosage \code{j} for SNP 2.
//' @param lg A vector with the log-probabilities of the genotypes for the first SNP. Element
//'     \code{j} is the log-probability of dosage \code{j - 1} in SNP 1.
//' @param lh A vector with the log-probabilities of the genotypes for the second SNP. Element
//'     \code{j} is the log-probability of dosage \code{j - 1} in SNP 2.
//'
//' @author David Gerard
//'
//' @seealso
//' \describe{
//'   \item{\code{\link{correst}}}{For the function that optimizes \code{corrlike_r}.}
//' }
//'
// [[Rcpp::export]]
double corrlike(double atanh_rho,
                const NumericMatrix& lX,
                const NumericMatrix& lY,
                const NumericVector& lg,
                const NumericVector& lh) {
  // check input -------------------------------------------
  int K = lg.length();
  int n = lX.nrow();
  if ((n != lX.nrow()) | (n != lY.nrow())) {
    Rcpp::stop("corrlike: lX and lY should have the same number of rows.");
  }
  if ((K != lX.ncol()) |
      (K != lY.ncol()) |
      (K != lg.length()) |
      (K != lh.length())) {
    Rcpp::stop("corrlike: lX.col(), lY.col(), lg.length(), and lh.lenght() should be equal.");
  }

  double rho = std::tanh(atanh_rho);

  // Get cumulative density function ----------------------
  NumericVector g_cdf = Rcpp::exp(log_cum_sum_exp(lg));
  NumericVector h_cdf = Rcpp::exp(log_cum_sum_exp(lh));

  // Get cdf matrix ---------------------------------------
  Rcpp::NumericMatrix cdf_mat(K, K);
  for (int j = 0; j < K; j++) {
    for (int ell = 0; ell < K; ell++) {
      cdf_mat(j, ell) = pnormcop(g_cdf(j), h_cdf(ell), rho);
    }
  }

  // Get prior_mat from cdf_mat ---------------------------
  Rcpp::NumericMatrix prior_mat(K, K);
  double non_log_prior = 0;
  for (int j = 0; j < K; j++) {
    for (int ell = 0; ell < K; ell++) {
      if ((ell != 0) & (j != 0)) {
        non_log_prior = cdf_mat(j, ell) -
          cdf_mat(j, ell - 1) -
          cdf_mat(j - 1, ell) +
          cdf_mat(j - 1, ell - 1);
        if (non_log_prior < TOL) {
          prior_mat(j, ell) = R_NegInf;
        } else {
          prior_mat(j, ell) = std::log(non_log_prior);
        }
      } else if (j != 0) {
        non_log_prior = cdf_mat(j, ell) - cdf_mat(j - 1, ell);
        if (non_log_prior < TOL) {
          prior_mat(j, ell) = R_NegInf;
        } else {
          prior_mat(j, ell) = std::log(non_log_prior);
        }
      } else if (ell != 0) {
        non_log_prior = cdf_mat(j, ell) - cdf_mat(j, ell - 1);
        if (non_log_prior < TOL) {
          prior_mat(j, ell) = R_NegInf;
        } else {
          prior_mat(j, ell) = std::log(non_log_prior);
        }
      } else { // both are zero
        non_log_prior = cdf_mat(j, ell);
        if (non_log_prior < TOL) {
          prior_mat(j, ell) = R_NegInf;
        } else {
          prior_mat(j, ell) = std::log(non_log_prior);
        }
      }
    }
  }

  // Get log-likelihood -----------------------------------
  double llike = 0.0;
  double ind_cont = R_NegInf;
  for (int i = 0; i < n; i++) {
    ind_cont = R_NegInf;
    for (int j = 0; j < K; j++) {
      for (int ell = 0; ell < K; ell++) {
        ind_cont = updog::log_sum_exp_2(ind_cont, prior_mat(j, ell) + lX(i, j) + lY(i, ell));
      }
    }
    llike = llike + ind_cont;
  }
  if ((llike == R_NegInf) | (llike == R_PosInf) | (llike == R_NaN) | (llike == R_NaReal)) {
    // List errout = List::create(Named("atanh_rho") = atanh_rho,
    //                            Named("lX") = lX,
    //                            Named("lY") = lY,
    //                            Named("lg") = lg,
    //                            Named("lh") = lh);
    Rcpp::Rcout << atanh_rho << std::endl;
    Rcpp::stop("corrlike: non-finite likelihood");
  }
  return llike;
}
