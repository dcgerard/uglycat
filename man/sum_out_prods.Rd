% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{sum_out_prods}
\alias{sum_out_prods}
\title{Calculate sum_n sum_j sum_l cross_mat(j, l) * pm1(n, j) * pm2(n, l)}
\usage{
sum_out_prods(pm1, pm2, crossmat)
}
\arguments{
\item{pm1}{A matrix of posterior probabilities for SNP 1. Element (i, j)
is the posterior probability of genotype j - 1 in individual i.}

\item{pm2}{A matrix of posterior probabilities for SNP 2. Element (i, j)
is the posterior probability of genotype j - 1 in individual i.}

\item{crossmat}{A matrix of probabilities. Element (j,k) is the
prior probability of genotype j-1 in SNP 1 and k-1 in SNP 2.}
}
\description{
Calculate sum_n sum_j sum_l cross_mat(j, l) * pm1(n, j) * pm2(n, l)
}
\author{
David Gerard
}
