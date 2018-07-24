context("log_sum_exp")

test_that("exp of log_cum_sum_exp is same as cumsum" {
  set.seed(1)
  x <- runif(11)
  expect_equal(cumsum(x), exp(log_cum_sum_exp(log(x))))
})

test_that("exp of log_sum_exp is same as sum", {
  set.seed(1)
  x <- runif(11)
  expect_equal(sum(x), exp(log_sum_exp(log(x))))
  expect_equal(log_sum_exp(c(-Inf, -Inf, -Inf)), -Inf)
})

test_that("exp of log_sum_exp_2 is same as sum", {
  set.seed(1)
  expect_equal(exp(log_sum_exp_2(log(2), log(4))), 6)
  expect_equal(log_sum_exp_2(-Inf, -Inf), -Inf)
})
