// Calculate genotype likelihoods

#include <RcppArmadillo.h>
#include <updog.h>
using namespace Rcpp;

//' Calculate the genotypes log-likelihoods.
//'
//' @param refvec A vector of reference counts.
//' @param sizevec A vector of size counts.
//' @param ploidy The ploidy of the species.
//' @param seq The sequencing error rate.
//' @param bias The allele bias parameter.
//' @param od The overdispersion parameter.
//'
//'
//' @return A matrix of numerics. The rows index the individuals and the
//'     columns index the genotype. Element (i, j) is the log-likelihood
//'     of dosage j for individual i.
//'
//' @author David Gerard
//'
//'
// [[Rcpp::export]]
NumericMatrix gen_like(NumericVector refvec,
                       NumericVector sizevec,
                       int ploidy,
                       double seq,
                       double bias,
                       double od) {
  // Check input -----------------------------------------------------------
  int nind = refvec.length();
  if (nind != sizevec.length()) {
    Rcpp::stop("get_wik_mat: sizevec and refvec must have the same length.");
  }

  // Calculate the posterior probability of each genotype -------------------
  NumericMatrix gen_like_mat(nind, ploidy + 1);
  NumericVector xi(ploidy + 1);
  for (int k = 0; k <= ploidy; k++) {
    xi(k) = updog::xi_double((double)k / (double)ploidy, seq, bias);
  }

  for (int i = 0; i < nind; i++) {
    for (int k = 0; k <= ploidy; k++) {
      gen_like_mat(i, k) = updog::dbetabinom_double(refvec(i), sizevec(i), xi(k), od, true);
    }
  }

  return gen_like_mat;
}
