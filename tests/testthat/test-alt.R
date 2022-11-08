test_that("alt gives same as hwep", {
  set.seed(1)
  x <- c(1, 2, 7, 5, 2)
  alpha <- rep(1, 5)
  ## exact
  t1 <- hwep:::ddirmult(x = x, alpha = alpha, lg = TRUE)
  ## approx
  t2 <- marg_alt_g(x = x, alpha = alpha, lg = TRUE, chains = 1)

  expect_equal(t1, t2, tolerance = 0.01)
})
