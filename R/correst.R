# Correlation estimation functions.

#' Estimate the correlation between two SNPs from the genotype likelihoods and the
#' genotype distribution.
#'
#' @param X A matrix of genotype likelihoods. Element \code{(i, j)} is the
#'     probability of the data for individual \code{i} given dosage \code{j - 1} for
#'     SNP 1.
#' @param Y A matrix of genotype likelihoods. Element \code{(i, j)} is the
#'     probability of the data for individual \code{i} given dosage \code{j - 1} for
#'     SNP 2.
#' @param g A vector of probabilities. Element \code{j} is the probability
#'     of dosage \code{j - 1} for SNP 1.
#' @param h A vector of probabilities. Element \code{j} is the probability
#'     of dosage \code{j - 1} for SNP 2.
#' @param is_log A logical. Are \code{X}, \code{Y}, \code{g}, and \code{h} all on the
#'     log-scale (\code{TRUE}) or not (\code{FALSE}).
#' @param use_cpp A logical Should we use the C++ implementation of the
#'     optimization (\code{TRUE}) or not (\code{FALSE})? The fully C++ implementation
#'     seems to be about 20 percent faster.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{zhat}}{The estimate of the atanh of the correlation of the SNPs on the
#'       transformed scale}
#'   \item{\code{se_zhat}}{The estimated standard error of zhat.}
#'   \item{\code{corest}}{The estimated correlation of the SNPs on the original scale.}
#' }
#'
#' @author David Gerard
#'
#' @export
#'
correst <- function(X, Y, g, h, is_log = FALSE, use_cpp = TRUE) {
  ## Test input ---------------------------------------
  assertthat::assert_that(is.matrix(X))
  assertthat::assert_that(is.numeric(X))
  assertthat::assert_that(is.matrix(Y))
  assertthat::assert_that(is.numeric(Y))
  assertthat::assert_that(is.vector(g, mode = "numeric"))
  assertthat::assert_that(is.vector(h, mode = "numeric"))
  assertthat::are_equal(nrow(X), nrow(Y))
  assertthat::are_equal(ncol(X), ncol(Y), length(g), length(h))
  assertthat::assert_that(is.vector(is_log, mode = "logical"))
  assertthat::are_equal(1, length(is_log))
  TOL <- 10 ^ -6
  if (!is_log) {
    stopifnot(abs(sum(g) - 1) < TOL, abs(sum(h) - 1) < TOL)
  } else {
    stopifnot(abs(updog::log_sum_exp(g)) < TOL, abs(updog::log_sum_exp(h)) < TOL)
  }

  ## Get log components -------------------------------
  if (is_log) {
    lX <- X
    lY <- Y
    lh <- h
    lg <- g
  } else {
    lX <- log(X)
    lY <- log(Y)
    lg <- log(g)
    lh <- log(h)
  }

  ## Initialize correlation via grid search
  nrho <- 5
  atanh_rhovec <- atanh(seq(-0.99, 0.99, length = nrho))
  lvec         <- rep(NA, length = length(atanh_rhovec))
  for (index in 1:length(atanh_rhovec)) {
    lvec[index] <- corrlike(atanh_rho = atanh_rhovec[index], lX = lX, lY = lY, lg = lg, lh = lh)
  }
  atanh_rho_init <- atanh_rhovec[which.max(lvec)]

  if (use_cpp) {
    ## Use new corr_optim function to get maximizer of corrlike ---
    oout <- corr_optim(atanh_rho = atanh_rho_init,
                       lX        = lX,
                       lY        = lY,
                       lg        = lg,
                       lh        = lh)
  } else{
    ## Use optim to get maximizer of corrlike -----------
    oout <- stats::optim(par     = atanh_rho_init,
                         fn      = corrlike,
                         gr      = NULL,
                         method  = "L-BFGS-B",
                         lower   = -6,
                         upper   = 6,
                         control = list(fnscale = -1),
                         hessian = TRUE,
                         lX      = lX,
                         lY      = lY,
                         lg      = lg,
                         lh      = lh)
  }

  ## Get parameter estimate and se -------------------------
  zhat <- oout$par
  se_zhat <- 1 / sqrt(max(-oout$hessian, 0))

  ## Get MLE of correlation of genotypes (not transformed genotypes)--------------------------------
  jd <- exp(dist_from_marg(g_cdf = exp(log_cum_sum_exp(lg)),
                           h_cdf = exp(log_cum_sum_exp(lh)),
                           rho = tanh(zhat)))
  sumjd <- sum(jd)
  jd <- jd / sumjd
  corest <- updog::oracle_cor_from_joint(jd)

  return(list(zhat = zhat, se_zhat = se_zhat, corest = corest))
}


#' Estimate correlation between two SNPs from two updog fits.
#'
#' @param uout1 A member of class \code{\link[updog]{flexdog}},
#'     which is the output of fitting \code{flexdog} to one SNP.
#' @param uout2 A member of class \code{\link[updog]{flexdog}},
#'     which is the output of fitting \code{flexdog} to the other SNP.
#' @param use_cpp A logical. Should we use the C++ implementation (\code{TRUE})
#'     or not (\code{FALSE})?
#'
#' @export
#'
#' @author David Gerard
correst_updog <- function(uout1, uout2, use_cpp = TRUE) {
  assertthat::assert_that(updog::is.flexdog(uout1))
  assertthat::assert_that(updog::is.flexdog(uout2))

  refvec1       <- uout1$input$refvec
  sizevec1      <- uout1$input$sizevec
  na1           <- is.na(refvec1) | is.na(sizevec1)
  refvec1[na1]  <- 0
  sizevec1[na1] <- 0

  refvec2       <- uout2$input$refvec
  sizevec2      <- uout2$input$sizevec
  na2           <- is.na(refvec2) | is.na(sizevec2)
  refvec2[na2]  <- 0
  sizevec2[na2] <- 0

  lX <- gen_like(refvec  = refvec1,
                 sizevec = sizevec1,
                 ploidy  = uout1$input$ploidy,
                 seq     = uout1$seq,
                 bias    = uout1$bias,
                 od      = uout1$od)
  lY <- gen_like(refvec  = refvec2,
                 sizevec = sizevec2,
                 ploidy  = uout2$input$ploidy,
                 seq     = uout2$seq,
                 bias    = uout2$bias,
                 od      = uout2$od)
  lg <- log(uout1$gene_dist)
  lh <- log(uout2$gene_dist)

  cout <- correst(X = lX, Y = lY, g = lg, h = lh, is_log = TRUE, use_cpp = use_cpp)

  return(cout)
}
