// Functions for use in correlation estimation

#include <Rcpp.h>
#include <math.h> // for M_PI
using namespace Rcpp;

// External functions ----------------------
double log_sum_exp(NumericVector x); // From utils.cpp
NumericVector log_cum_sum_exp(NumericVector x); // From utils.cpp
double log_sum_exp_2(double x, double y); // From utils.cpp
double expit(double x); // From utils.cpp
double pnormcop(double x, double y, double rho); // From gauss_cop.cpp

// Global variables -------------------------
double TOL = 1000.0 * DBL_EPSILON;

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
  for (int j = 0; j < K; j++) {
    for (int ell = 0; ell < K; ell++) {
      if ((ell != 0) & (j != 0)) {
        prior_mat(j, ell) = std::log(cdf_mat(j, ell) -
          cdf_mat(j, ell - 1) -
          cdf_mat(j - 1, ell) +
          cdf_mat(j - 1, ell - 1));
      } else if (j != 0) {
        prior_mat(j, ell) = std::log(cdf_mat(j, ell) - cdf_mat(j - 1, ell));
      } else if (ell != 0) {
        prior_mat(j, ell) = std::log(cdf_mat(j, ell) - cdf_mat(j, ell - 1));
      } else { // both are zero
        prior_mat(j, ell) = std::log(cdf_mat(j, ell));
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
        ind_cont = log_sum_exp_2(ind_cont, prior_mat(j, ell) + lX(i, j) + lY(i, ell));
      }
    }
    llike = llike + ind_cont;
  }
  return llike;
}
