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
#' @param method Which method should we use to estimate the correlation? MLE using
#'     C++ (\code{"mleCpp"}) or  MLE using R (\code{"mleR"})?
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{zhat_gc}}{The estimate of the atanh of the correlation of the SNPs on the
#'       gaussian copula scale}
#'   \item{\code{sezhat_gc}}{The estimated standard error of zhat_gc.}
#'   \item{\code{corest}}{The estimated correlation of the SNPs on the original scale.}
#'   \item{\code{zhat}}{\code{atanh(corest)}}
#'   \item{\code{sezhat}}{The estimated standard error of zhat.}
#' }
#'
#' @author David Gerard
#'
#' @export
#'
correst <- function(X,
                    Y,
                    g,
                    h,
                    is_log = FALSE,
                    method = c("mleCpp", "mleR")) {
  ## Test input ---------------------------------------
  method <- match.arg(method)
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

  ## Check for monomorphic snps --------------------------------------------
  if ((abs(max(exp(lg)) - 1) < TOL) | (abs(max(exp(lh)) - 1) < TOL)) {
    return(list(zhat_gc   = NA,
                sezhat_gc = NA,
                corest    = NA,
                zhat      = NA,
                sezhat    = NA,
                comment   = "One SNP is monomorphic."))
  }

  ## Initialize correlation via grid search --------------------------------
  nrho <- 5
  atanh_rhovec <- atanh(seq(-0.99, 0.99, length = nrho))
  lvec         <- rep(NA, length = length(atanh_rhovec))
  for (index in 1:length(atanh_rhovec)) {
    lvec[index] <- corrlike(atanh_rho = atanh_rhovec[index], lX = lX, lY = lY, lg = lg, lh = lh)
  }
  atanh_rho_init <- atanh_rhovec[which.max(lvec)]

  if (method == "mleCpp") {
    ## Use new corr_optim function to get maximizer of corrlike ---------
    oout <- corr_optim(atanh_rho = atanh_rho_init,
                       lX        = lX,
                       lY        = lY,
                       lg        = lg,
                       lh        = lh)
  } else if (method == "mleR") {
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
  } else {
    stop("correst: how did you get here?")
  }

  ## Get parameter estimate and se on gaussian copula scale -------------------------
  zhat <- oout$par
  sezhat_gc <- 1 / sqrt(max(-oout$hessian, 0))

  ## Get MLE of correlation of genotypes (not transformed genotypes)--------------------------------
  g_cdf <- exp(log_cum_sum_exp(lg))
  h_cdf <- exp(log_cum_sum_exp(lh))
  corest <- z_to_corr(z = zhat, g_cdf = g_cdf, h_cdf = h_cdf)

  ## Get zscores and ses on original scale -------------------------------
  myenv <- new.env()
  assign(x = "zhat", value = zhat, envir = myenv)
  assign(x = "g_cdf", value = g_cdf, envir = myenv)
  assign(x = "h_cdf", value = h_cdf, envir = myenv)
  nout <- stats::numericDeriv(quote(z_to_z(z = zhat, g_cdf = g_cdf, h_cdf = h_cdf)), "zhat", myenv)
  sezhat <- sezhat_gc * abs(c(attr(nout, "gradient")))

  return(list(zhat_gc = zhat,
              se_zhat = sezhat_gc,
              corest  = corest,
              zhat    = atanh(corest),
              sezhat  = sezhat))
}


#' Estimate correlation between two SNPs from two updog fits.
#'
#' @param uout1 A member of class \code{\link[updog]{flexdog}},
#'     which is the output of fitting \code{flexdog} to one SNP.
#' @param uout2 A member of class \code{\link[updog]{flexdog}},
#'     which is the output of fitting \code{flexdog} to the other SNP.
#' @param method Which method should we use to estimate the correlation? MLE using
#'     C++ (\code{"mleCpp"}) or  MLE using R (\code{"mleR"})?
#'
#' @export
#'
#' @author David Gerard
correst_updog <- function(uout1, uout2, method = c("mleCpp", "mleR")) {
  method <- match.arg(method)
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

  cout <- correst(X = lX, Y = lY, g = lg, h = lh, is_log = TRUE, method = method)

  return(cout)
}


#' Estimate correlation using weights that assume independence.
#'
#' @param pm1 The posterior probability matrix for SNP 1. Element (i, j) is
#'     the posterior probability that individual i has genotype j - 1.
#' @param pm2 The posterior probability matrix for SNP 2. Element (i, j) is
#'     the posterior probability that individual i has genotype j - 1.
#' @param g The pmf of SNP 1. ELement j is the posterior probability of genotype j - 1.
#' @param h The pmf of SNP 2. ELement j is the posterior probability of genotype j - 1.
#'
#'
#' @author David Gerard
correst_ind <- function(pm1, pm2, g, h) {
  ## Check input --------------------------------------
  assertthat::assert_that(is.matrix(pm1))
  assertthat::assert_that(is.numeric(pm1))
  assertthat::assert_that(is.matrix(pm2))
  assertthat::assert_that(is.numeric(pm2))
  assertthat::assert_that(is.vector(g, mode = "numeric"))
  assertthat::assert_that(is.vector(h, mode = "numeric"))
  assertthat::are_equal(length(g), length(h), ncol(pm1), ncol(pm2))
  assertthat::are_equal(nrow(pm1), nrow(pm2))

  ## Get correlation estimate ------------------
  ploidy <- length(g) - 1
  assertthat::assert_that(ploidy >= 1)

  mug   <- sum(0:ploidy * g)
  sig2g <- sum((0:ploidy - mug) ^ 2 * g)
  muh   <- sum(0:ploidy * h)
  sig2h <- sum((0:ploidy - muh) ^ 2 * h)

  crossmat <- outer(0:ploidy - mug, 0:ploidy - muh, "*")
  navec <- !(is.na(pm1[, 1]) | is.na(pm2[, 1]))

  covout = sum_out_prods(pm1 = pm1[navec, ], pm2 = pm2[navec, ], crossmat = crossmat)

  corrout <- covout / (nrow(pm1) * sqrt(sig2g) * sqrt(sig2h))

  return(corrout)
}


#' Estimate correlation assuming independent weights using updog output.
#'
#' @inheritParams correst_updog
#'
#' @author David Gerard
#'
correst_ind_updog <- function(uout1, uout2) {
  assertthat::assert_that(updog::is.flexdog(uout1))
  assertthat::assert_that(updog::is.flexdog(uout2))

  corrout <- correst_ind(pm1 = uout1$postmat,
                         pm2 = uout2$postmat,
                         g = uout1$gene_dist,
                         h = uout2$gene_dist)

  return(corrout)
}


#' Transforms from z of gaussian copula to correlation of multivariate distribution.
#'
#' @param z Any real number. The atanh of the underlying gaussian
#'     copula correlation.
#' @param g_cdf One marginal cdf. Element i is the probability of
#'    of being less than or equal to i - 1.
#' @param h_cdf The other marginal cdf.Element i is the probability of
#'    of being less than or equal to i - 1.
#'
#' @author David Gerard
z_to_corr <- function(z, g_cdf, h_cdf) {
  jd <- exp(dist_from_marg(g_cdf = g_cdf,
                           h_cdf = h_cdf,
                           rho = tanh(z)))
  sumjd <- sum(jd)
  jd <- jd / sumjd
  correst <- updog::oracle_cor_from_joint(jd)
  return(correst)
}

z_to_z <- function(z, g_cdf, h_cdf) {
  atanh(z_to_corr(z, g_cdf, h_cdf))
}

#' Estimate correlation using one Newton step from result of correst_ind_updog
#'
#' @inheritParams correst_updog
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{zhat_gc}}{The estimate of the atanh of the correlation of
#'       the SNPs on the gaussian copula scale}
#'   \item{\code{sezhat_gc}}{The estimated standard error of zhat_gc.}
#'   \item{\code{corest}}{The estimated correlation of the SNPs on the
#'       original scale.}
#'   \item{\code{zhat}}{\code{atanh(corest)}}
#'   \item{\code{sezhat}}{The estimated standard error of zhat.}
#' }
#'
#' @author David Gerard
#'
#' @export
#'
correst_onestep_updog <- function(uout1, uout2) {
  assertthat::assert_that(updog::is.flexdog(uout1))
  assertthat::assert_that(updog::is.flexdog(uout2))

  ## Check to make sure neither is degenerate --------
  TOL <- 10 ^ -6
  if ((abs(max(uout1$gene_dist) - 1) < TOL) | (abs(max(uout2$gene_dist) - 1) < TOL)) {
    return(list(zhat_gc   = NA,
                sezhat_gc = NA,
                corest    = NA,
                zhat      = NA,
                sezhat    = NA,
                comment   = "One SNP is monomorphic."))
  }

  ## Initial naive corrlation estimate --------------
  corrhat_naive <- correst_ind_updog(uout1 = uout1, uout2 = uout2)

  ## Transform to Gaussian copula space -------------
  g_cdf <- cumsum(uout1$gene_dist)
  h_cdf <- cumsum(uout2$gene_dist)
  diffz <- function(z) {
    corrhat_naive - z_to_corr(z = z, g_cdf = g_cdf, h_cdf = h_cdf)
  }
  if (sign(diffz(-6) != sign(diffz(6)))) {
    uniout <- stats::uniroot(f = diffz, interval = c(-6, 6))
    ztilde <- uniout$root
  } else {
    ztilde <- atanh(corrhat_naive)
  }


  ## Get genotype likelihoods ------------------------
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

  ## One Newton step to get final est --------------
  eps <- 0.001
  like_meps <- corrlike(atanh_rho = ztilde - eps, lX = lX, lY = lY, lg = lg, lh = lh)
  like_peps <- corrlike(atanh_rho = ztilde + eps, lX = lX, lY = lY, lg = lg, lh = lh)
  like_0 <- corrlike(atanh_rho = ztilde, lX = lX, lY = lY, lg = lg, lh = lh)

  diff1 <- (like_peps - like_meps) / (2 * eps)
  diff2 <- (like_peps - 2 * like_0 + like_meps) / (eps^2)
  zfinal <- ztilde - diff1 / diff2

  ## Get Hessian again ------------------------------------------
  like_meps <- corrlike(atanh_rho = zfinal - eps, lX = lX, lY = lY, lg = lg, lh = lh)
  like_peps <- corrlike(atanh_rho = zfinal + eps, lX = lX, lY = lY, lg = lg, lh = lh)
  like_0 <- corrlike(atanh_rho = zfinal, lX = lX, lY = lY, lg = lg, lh = lh)
  diff2 <- (like_peps - 2 * like_0 + like_meps) / (eps^2)

  sezhat_gc <- 1 / sqrt(max(-diff2, 0))

  ## Get final correlation --------------------------------------
  corest <- z_to_corr(z = zfinal, g_cdf = g_cdf, h_cdf = h_cdf)

  ## Get zscores and ses on original scale -------------------------------
  myenv <- new.env()
  assign(x = "zfinal", value = zfinal, envir = myenv)
  assign(x = "g_cdf", value = g_cdf, envir = myenv)
  assign(x = "h_cdf", value = h_cdf, envir = myenv)
  nout <- stats::numericDeriv(quote(z_to_z(z = zfinal, g_cdf = g_cdf, h_cdf = h_cdf)), "zfinal", myenv)
  sezhat <- sezhat_gc * abs(c(attr(nout, "gradient")))

  return(list(zhat_gc   = zfinal,
              sezhat_gc = sezhat_gc,
              corest    = corest,
              zhat      = atanh(corest),
              sezhat    = sezhat))
}

