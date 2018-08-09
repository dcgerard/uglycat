context("corrtest")

test_that("corrtest works", {
  set.seed(1)
  itermax <- 500
  n <- 100
  k <- 7
  X <- matrix(stats::runif(n = n * k), nrow = n, ncol = k)
  Y <- matrix(stats::runif(n = n * k), nrow = n, ncol = k)
  g <- stats::runif(k)
  g <- g / sum(g)
  h <- g
  ceout <- correst(X = X, Y = Y, g = g, h = h)

  rhovec <- seq(-0.99, 0.99, length = 100)
  lvec <- rep(NA, length = length(rhovec))
  for (index in 1:length(rhovec)) {
    lvec[index] <- corrlike(atanh_rho = atanh(rhovec[index]), lX = log(X), lY = log(Y), lg = log(g), lh = log(h))
  }
})

test_that("corroptim works", {
  if (.Platform$OS.type == "windows") {
    skip("optim is throwing weird results in i386.")
  }
  set.seed(1)
  itermax <- 500
  n <- 100
  k <- 7
  X <- matrix(stats::runif(n = n * k), nrow = n, ncol = k)
  Y <- matrix(stats::runif(n = n * k), nrow = n, ncol = k)
  g <- stats::runif(k)
  g <- g / sum(g)
  h <- g

  lX <- log(X)
  lY <- log(Y)
  lg <- log(g)
  lh <- log(h)

  atanh_rho <- 5

  coout <- corr_optim(atanh_rho = atanh_rho,
                      lX = lX,
                      lY = lY,
                      lg = lg,
                      lh = lh)

  oout <- optim(par = atanh_rho,
                fn = corrlike,
                method = "L-BFGS-B",
                lower = -6,
                upper = 6,
                control = list(fnscale = -1),
                hessian = TRUE,
                lX = lX,
                lY = lY,
                lg = lg,
                lh = lh)

  expect_equal(oout$par, coout$par, tol = 10^-3)
  expect_equal(oout$value, coout$value, tol = 10 ^ -3)
  expect_equal(c(oout$hessian), coout$hessian, tol = 10^-3)

  # two times speedup
  # microbenchmark::microbenchmark(
  #   coout <- corr_optim(atanh_rho = atanh_rho,
  #                       lX = lX,
  #                       lY = lY,
  #                       lg = lg,
  #                       lh = lh),
  #   oout <- optim(par = atanh_rho,
  #                 fn = corrlike,
  #                 method = "L-BFGS-B",
  #                 lower = -6,
  #                 upper = 6,
  #                 control = list(fnscale = -1),
  #                 hessian = TRUE,
  #                 lX = lX,
  #                 lY = lY,
  #                 lg = lg,
  #                 lh = lh)
  # )

})

test_that("correst works on real data", {
  udat <- readRDS("udat.RDS")
  cout <- correst_updog(uout1 = udat$uout1, uout2 = udat$uout2)
  coutind <- correst_ind_updog(uout1 = udat$uout1, uout2 = udat$uout2)
  coutone <- correst_onestep_updog(uout1 = udat$uout1, uout2 = udat$uout2)
})

test_that("correst works on real data 2", {
  udat <- readRDS("harddat.RDS")
  cout <- correst_updog(uout1 = udat$uout1, uout2 = udat$uout2)
  coutind <- correst_ind_updog(uout1 = udat$uout1, uout2 = udat$uout2)
  coutone <- correst_onestep_updog(uout1 = udat$uout1, uout2 = udat$uout2)

  # microbenchmark::microbenchmark(
  #   cout <- correst_updog(uout1 = udat$uout1, uout2 = udat$uout2, method = "mleR"),
  #   cout <- correst_updog(uout1 = udat$uout1, uout2 = udat$uout2, method = "mleCpp"),
  #   cout <- correst_ind_updog(uout1 = udat$uout1, uout2 = udat$uout2),
  #   cout <- correst_onestep_updog(uout1 = udat$uout1, uout2 = udat$uout2)
  # )
})

test_that("correst works on real data 3", {
  udat <- readRDS("harddat2.RDS")
  cout <- correst_updog(uout1 = udat$uout1, uout2 = udat$uout2)
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

test_that("errdat.RDS does not return NaN", {
  errdat <- readRDS("./errdat.RDS")
  clout <- corrlike(atanh_rho = errdat$atanh_rho,
                    lX = errdat$lX,
                    lY = errdat$lY,
                    lg = errdat$lg,
                    lh = errdat$lh)
  expect_true(!is.nan(clout))
})
