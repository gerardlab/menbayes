test_that("alt gives same as hwep", {
  set.seed(1)
  x <- c(1, 2, 7, 5, 2)
  alpha <- rep(1, 5)
  ## exact
  t1 <- hwep:::ddirmult(x = x, alpha = alpha, lg = TRUE)
  ## approx
  t2 <- marg_alt_g(x = x, beta = alpha, lg = TRUE, chains = 1)

  expect_equal(t1, t2, tolerance = 0.01)
})

test_that("alt_gl gives same as alt_g when have perfect gl", {
  set.seed(1)
  x <- c(1, 2, 4, 2, 1)
  ploidy <- length(x) - 1
  nind <- sum(x)
  alpha <- rep(1, 5)

  t1 <- marg_alt_g(x = x, beta = alpha, lg = TRUE, chains = 1)

  gl <- matrix(-Inf, nrow = nind, ncol = ploidy + 1)
  genovec <- unlist(mapply(0:ploidy, x, FUN = rep))

  for (i in seq_len(nind)) {
    gl[i, genovec[[i]] + 1] <- 0
  }

  t2 <- marg_alt_gl(gl = gl, beta = alpha, lg = TRUE, chains = 1)

  ## Should differ by multinomial normalizing constant since in t1
  ## individuals are indistinguishable, but in t2 they are.
  normconst <- lfactorial(nind) - sum(lfactorial(x))

  expect_equal(t1, t2 + normconst, tolerance = 0.1)
})
