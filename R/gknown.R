###################
## Methods when genotypes are known
###################

#' Stan version of marg_alt_g(). Not to be used.
#'
#' @noRd
marg_alt_g_stan <- function(x, beta = rep(1, length(x)), lg = TRUE, ...) {
  ploidy <- length(x) - 1
  stopifnot(length(x) == length(beta))
  stan_dat <- list(K = ploidy, x = x, beta = beta)
  stan_out <- rstan::sampling(object = stanmodels$alt_g,
                              data = stan_dat,
                              verbose = FALSE,
                              show_messages = FALSE,
                              ...)
  bridge_out <- bridgesampling::bridge_sampler(stan_out, verbose = FALSE, silent = TRUE)

  if (lg) {
    mx <- bridge_out$logml
  } else {
    mx <- exp(bridge_out$logml)
  }
  return(mx)
}

#' Density of dirichlet-multinomial
#'
#' @param x the counts
#' @param alpha the dirichlet parameters
#' @param lg should we log marginal likelihood?
#'
#' @noRd
ddirmult <- function (x, alpha, lg = FALSE)
{
    stopifnot(length(x) == length(alpha))
    stopifnot(all(alpha > 0))
    asum <- sum(alpha)
    n <- sum(x)
    ll <- lgamma(asum) +
      lgamma(n + 1) -
      lgamma(n + asum) +
      sum(lgamma(x + alpha)) -
      sum(lgamma(alpha)) -
      sum(lgamma(x + 1))
    if (!lg) {
        ll <- exp(ll)
    }
    return(ll)
}

#' Marginal likelihood under alternative when genotypes are known
#'
#' @param x The genotype counts
#' @param beta The prior hyperparameters
#' @param lg A logical. Should we log the marginal likelihood or not?
#' @param ... Additional arguments to pass to \code{\link[rstan]{sampling}()}.
#'
#' @return The marginal likelihood.
#'
#' @author Mira Thakkar and David Gerard
#'
#' @examples
#' x <- c(3L, 21L, 52L, 20L, 4L)
#' marg_alt_g(x = x, chains = 1)
#'
#' @export
marg_alt_g <- function(x, beta = rep(1, length(x)), lg = TRUE, ...) {
  ploidy <- length(x) - 1
  stopifnot(length(x) == length(beta))
  mx <- ddirmult(x = x, alpha = beta, lg = lg)
  return(mx)
}

#' Marginal likelihood, no double reduction, no preferential pairing, genotypes known.
#'
#' @inheritParams marg_f1_dr_pp_g4
#'
#' @author Mira Thakkar and David Gerard
#'
#' @export
#'
#' @examples
#' x <- c(3L, 21L, 52L, 20L, 4L) ## simulated via null
#' g1 <- 2
#' g2 <- 2
#' mnull <- marg_f1_ndr_npp_g4(x = x, g1 = g1, g2 = g2)
#' malt <- -15.34 ## output of marg_alt_g(x = x)
#' mnull - malt ## log-BF
#'
#' x <- c(25L, 24L, 20L, 15L, 16L) ## simulated via alt
#' g1 <- 2
#' g2 <- 2
#' mnull <- marg_f1_ndr_npp_g4(x = x, g1 = g1, g2 = g2)
#' malt <- -15.34 ## output of marg_alt_g(x = x)
#' mnull - malt ## log-BF
#'
marg_f1_ndr_npp_g4 <- function(x,
                               g1,
                               g2,
                               lg = TRUE) {
  gf <- offspring_gf(alpha = 0, xi = 1/3, p1 = g1, p2 = g2)
  return(stats::dmultinom(x = x, prob = gf, log = lg))
}

#' Marginal likelihood, double reduction, no preferential pairing, genotypes known.
#'
#' Here, parental genotypes are known.
#'
#' @inheritParams marg_f1_dr_pp_g4
#'
#' @author Mira Thakkar and David Gerard
#'
#' @examples
#' x <- c(3L, 21L, 52L, 20L, 4L) ## simulated via null
#' g1 <- 2
#' g2 <- 2
#' mnull <- marg_f1_dr_npp_g4(x = x, g1 = g1, g2 = g2)
#' malt <- -15.34 ## output of marg_alt_g(x = x)
#' mnull - malt ## log-BF
#'
#' x <- c(25L, 24L, 20L, 15L, 16L) ## simulated via alt
#' g1 <- 2
#' g2 <- 2
#' mnull <- marg_f1_dr_npp_g4(x = x, g1 = g1, g2 = g2)
#' malt <- -15.34 ## output of marg_alt_g(x = x)
#' mnull - malt ## log-BF
#'
#' @export
marg_f1_dr_npp_g4 <- function(x,
                              g1,
                              g2,
                              lg = TRUE,
                              ...) {
  stopifnot(length(x) == 5,
            length(g1) == 1,
            length(g2) == 1)
  stopifnot(g1 >= 0, g1 <= 4, g2 >= 0, g2 <= 4)
  drbound <- hwep::drbounds(ploidy = 4)
  stan_dat <- list(x = x,
                   drbound = drbound,
                   g1 = g1,
                   g2 = g2)
  stan_out <- rstan::sampling(object = stanmodels$marg_dr_npp_g4,
                              data = stan_dat,
                              verbose = FALSE,
                              show_messages = FALSE,
                              ...)
  bridge_out <- bridgesampling::bridge_sampler(stan_out, verbose = FALSE, silent = TRUE)

  if (lg) {
    mx <- bridge_out$logml
  } else {
    mx <- exp(bridge_out$logml)
  }
  return(mx)
}

#' Marginal likelihood, no double reduction, preferential pairing, genotypes known.
#'
#' @inheritParams marg_f1_dr_pp_g4
#'
#' @author Mira Thakkar and David Gerard
#'
#' @examples
#' x <- c(3L, 21L, 52L, 20L, 4L) ## simulated via null
#' g1 <- 2
#' g2 <- 2
#' mnull <- marg_f1_ndr_pp_g4(x = x, g1 = g1, g2 = g2)
#' malt <- -15.34 ## output of marg_alt_g(x = x)
#' mnull - malt ## log-BF
#'
#' x <- c(25L, 24L, 20L, 15L, 16L) ## simulated via alt
#' g1 <- 2
#' g2 <- 2
#' mnull <- marg_f1_ndr_pp_g4(x = x, g1 = g1, g2 = g2)
#' malt <- -15.34 ## output of marg_alt_g(x = x)
#' mnull - malt ## log-BF
#'
#' @export
#'
marg_f1_ndr_pp_g4 <- function(x,
                              g1,
                              g2,
                              lg = TRUE,
                              ...) {
  stopifnot(length(x) == 5,
            length(g1) == 1,
            length(g2) == 1)
  stopifnot(g1 >= 0, g1 <= 4, g2 >= 0, g2 <= 4)
  stan_dat <- list(x = x,
                   g1 = g1,
                   g2 = g2)
  stan_out <- rstan::sampling(object = stanmodels$marg_ndr_pp_g4,
                              data = stan_dat,
                              verbose = FALSE,
                              show_messages = FALSE,
                              ...)
  bridge_out <- bridgesampling::bridge_sampler(stan_out, verbose = FALSE, silent = TRUE)

  if (lg) {
    mx <- bridge_out$logml
  } else {
    mx <- exp(bridge_out$logml)
  }
  return(mx)
}

#' Marginal likelihood in tetraploid F1 population
#'
#' Here, parental genotypes are known.
#'
#' @param x The genotype counts of the offspring. \code{x[i]} is the
#'     number of offspring with genotype \code{i-1}.
#' @param g1 The first parent's genotype.
#' @param g2 The second parent's genotype.
#' @param lg A logical. Should we log the marginal likelihood (\code{TRUE})
#'     or not (\code{FALSE})?
#' @param ... Additional arguments to pass to \code{\link[stan]{sampling}()}.
#'
#' @author Mira Thakkar and David Gerard
#'
#' @examples
#' x <- c(3L, 21L, 52L, 20L, 4L) ## simulated via null
#' g1 <- 2
#' g2 <- 2
#' mnull <- marg_f1_dr_pp_g4(x = x, g1 = g1, g2 = g2)
#' malt <- -15.34 ## output of marg_alt_g(x = x)
#' mnull - malt ## log-BF
#'
#' x <- c(25L, 24L, 20L, 15L, 16L) ## simulated via alt
#' g1 <- 2
#' g2 <- 2
#' mnull <- marg_f1_dr_pp_g4(x = x, g1 = g1, g2 = g2)
#' malt <- -15.34 ## output of marg_alt_g(x = x)
#' mnull - malt ## log-BF
#'
#' @export
marg_f1_dr_pp_g4 <- function(x,
                             g1,
                             g2,
                             lg = TRUE,
                             ...) {
  stopifnot(length(x) == 5,
            length(g1) == 1,
            length(g2) == 1)
  stopifnot(g1 >= 0, g1 <= 4, g2 >= 0, g2 <= 4)
  drbound <- hwep::drbounds(ploidy = 4)
  stan_dat <- list(x = x,
                   drbound = drbound,
                   g1 = g1,
                   g2 = g2)
  stan_out <- rstan::sampling(object = stanmodels$marg_dr_pp_g4,
                              data = stan_dat,
                              verbose = FALSE,
                              show_messages = FALSE,
                              ...)
  bridge_out <- bridgesampling::bridge_sampler(stan_out, verbose = FALSE, silent = TRUE)

  if (lg) {
    mx <- bridge_out$logml
  } else {
    mx <- exp(bridge_out$logml)
  }
  return(mx)
}
