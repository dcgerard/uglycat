context("corrtest")

test_that("corrtest works", {
  set.seed(1)
  itermax <- 500
  z <- rep(NA, length = itermax)
  sez <- rep(NA, length = itermax)
  # for (index in 1:itermax) {
  # n <- 100
  # k <- 7
  # X <- matrix(stats::runif(n = n * k), nrow = n, ncol = k)
  # Y <- matrix(stats::runif(n = n * k), nrow = n, ncol = k)
  # g <- stats::runif(k)
  # g <- g / sum(g)
  # h <- g
  # ceout <- correst(X = X, Y = Y, g = g, h = h)
  # z[index] <- ceout$zhat
  # sez[index] <- ceout$se_zhat
  # }
  # which_keep <- !is.na(z) & !is.na(sez)
  # z <- z[which_keep]
  # sez <- sez[which_keep]
  # library(ashr)
  # ash(betahat = z, sebetahat = sez) -> aout
  #
  #
  # rhovec <- seq(-0.99, 0.99, length = 100)
  # lvec <- rep(NA, length = length(rhovec))
  # for (index in 1:length(rhovec)) {
  #   lvec[index] <- corrlike(atanch_rho = atanh(rhovec[index]), lX = log(X), lY = log(Y), lg = log(g), lh = log(h))
  # }
  # plot(rhovec, lvec, type = "l")
})

test_that("correst works on real data", {
  library(updog)
  data(snpdat)
  snp1dat <- snpdat[snpdat$snp == "SNP2", ]
  snp2dat <- snpdat[snpdat$snp == "SNP3", ]
  uout1 <- flexdog(refvec  = snp1dat$counts,
                   sizevec = snp1dat$size,
                   ploidy  = 6,
                   model   = "s1",
                   bias    = 1,
                   verbose = FALSE)
  uout2 <- flexdog(refvec  = snp2dat$counts,
                   sizevec = snp2dat$size,
                   ploidy  = 6,
                   model   = "s1",
                   bias    = 1,
                   verbose = FALSE)
  cout <- correst_updog(uout1 = uout1, uout2 = uout2)
})

test_that("bvnl works", {
  if (requireNamespace("mvtnorm", quietly = TRUE)) {

    rho <- 0.5
    Rmat <- matrix(c(1, rho, rho, 1), ncol = 2, nrow = 2)
    mvval <- mvtnorm::pmvnorm(upper = c(2, 2), corr = Rmat)[[1]]
    uval <- bvnl(dh = 2, dk = 2, r = rho)
    expect_equal(mvval, uval, tol = 10^-5)

    rho <- -0.5
    Rmat <- matrix(c(1, rho, rho, 1), ncol = 2, nrow = 2)
    mvval <- mvtnorm::pmvnorm(upper = c(2, 2), corr = Rmat)[[1]]
    uval <- bvnl(dh = 2, dk = 2, r = rho)
    expect_equal(mvval, uval, tol = 10^-5)

    rho <- 0.926
    Rmat <- matrix(c(1, rho, rho, 1), ncol = 2, nrow = 2)
    mvval <- mvtnorm::pmvnorm(upper = c(2, 1), corr = Rmat)[[1]]
    uval <- bvnl(dh = 2, dk = 1, r = rho)
    cat(bvnl(dh = 2, dk = 2, r = 0.924), "\n")
    cat(bvnl(dh = 2, dk = 2, r = 0.926), "\n")
    expect_equal(mvval, uval, tol = 10^-5)

    rho <- -0.99
    Rmat <- matrix(c(1, rho, rho, 1), ncol = 2, nrow = 2)
    mvval <- mvtnorm::pmvnorm(upper = c(2, 2), corr = Rmat)[[1]]
    uval <- bvnl(dh = 2, dk = 2, r = rho)
    expect_equal(mvval, uval, tol = 10^-5)

    rho <- 1
    Rmat <- matrix(c(1, rho, rho, 1), ncol = 2, nrow = 2)
    mvval <- mvtnorm::pmvnorm(upper = c(2, 2), corr = Rmat)[[1]]
    uval <- bvnl(dh = 2, dk = 2, r = rho)
    expect_equal(mvval, uval, tol = 10^-5)

    rho <- -1
    Rmat <- matrix(c(1, rho, rho, 1), ncol = 2, nrow = 2)
    mvval <- mvtnorm::pmvnorm(upper = c(2, 2), corr = Rmat)[[1]]
    uval <- bvnl(dh = 2, dk = 2, r = rho)
    expect_equal(mvval, uval, tol = 10^-5)

    rho <- 0
    Rmat <- matrix(c(1, rho, rho, 1), ncol = 2, nrow = 2)
    mvval <- mvtnorm::pmvnorm(upper = c(2, 2), corr = Rmat)[[1]]
    uval <- bvnl(dh = 2, dk = 2, r = rho)
    expect_equal(mvval, uval, tol = 10^-5)

    rho <- 0.1
    Rmat <- matrix(c(1, rho, rho, 1), ncol = 2, nrow = 2)
    mvval <- mvtnorm::pmvnorm(upper = c(2, 2), corr = Rmat)[[1]]
    uval <- bvnl(dh = 2, dk = 2, r = rho)
    expect_equal(mvval, uval, tol = 10^-5)

    rho <- -0.1
    Rmat <- matrix(c(1, rho, rho, 1), ncol = 2, nrow = 2)
    mvval <- mvtnorm::pmvnorm(upper = c(2, 2), corr = Rmat)[[1]]
    uval <- bvnl(dh = 2, dk = 2, r = rho)
    expect_equal(mvval, uval, tol = 10^-5)

    rho <- 0.5
    Rmat <- matrix(c(1, rho, rho, 1), ncol = 2, nrow = 2)
    mvval <- mvtnorm::pmvnorm(upper = c(2, Inf), corr = Rmat)[[1]]
    uval <- bvnl(dh = 2, dk = Inf, r = rho)
    expect_equal(mvval, uval, tol = 10^-5)

    rho <- 0.5
    Rmat <- matrix(c(1, rho, rho, 1), ncol = 2, nrow = 2)
    mvval <- mvtnorm::pmvnorm(upper = c(-Inf, Inf), corr = Rmat)[[1]]
    uval <- bvnl(dh = -Inf, dk = Inf, r = rho)
    expect_equal(mvval, uval, tol = 10^-5)

    rho <- 0.5
    Rmat <- matrix(c(1, rho, rho, 1), ncol = 2, nrow = 2)
    mvval <- mvtnorm::pmvnorm(upper = c(Inf, Inf), corr = Rmat)[[1]]
    uval <- bvnl(dh = Inf, dk = Inf, r = rho)
    expect_equal(mvval, uval, tol = 10^-5)

    rho <- 0.5
    Rmat <- matrix(c(1, rho, rho, 1), ncol = 2, nrow = 2)
    mvval <- mvtnorm::pmvnorm(upper = c(-Inf, -Inf), corr = Rmat)[[1]]
    uval <- bvnl(dh = -Inf, dk = -Inf, r = rho)
    expect_equal(mvval, uval, tol = 10^-5)

    # rho <- 0.5
    # Rmat <- matrix(c(1, rho, rho, 1), ncol = 2, nrow = 2)
    # microbenchmark::microbenchmark(mvval <- mvtnorm::pmvnorm(upper = c(2, 2), corr = Rmat)[[1]],
    #                                uval <- bvnl(dh = 2, dk = 2, r = rho))
  }
})

test_that("pnorm_cop_works", {
  expect_equal(pnormcop(x = 0.5, y = 1, rho = 0.4), 0.5)
  expect_equal(pnormcop(x = 1, y = 0.5, rho = 0.4), 0.5)
  expect_equal(pnormcop(x = 0.2, y = 0.2, rho = 0), 0.04)
  expect_equal(pnormcop(x = 0, y = 1, rho = 0.4), 0)
})
