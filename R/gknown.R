###################
## Methods when genotypes are known
###################

#' Marginal likelihood under alternative when genotypes are known
#'
#' @param x The genotype counts
#' @param beta The prior hyperparameters
#' @param lg A logical. Should we log the marginal likelihood or not?
#' @param ... Additional arguments to pass to \code{\link[rstan]{sampling}()}.
#'
#' @return The marginal likelihood.
#'
#' @author David Gerard
#'
#' @examples
#' x <- c(1, 2, 7, 3, 3)
#' beta <- rep(1, 5)
#' marg_alt_g(x = x, beta = beta, chains = 1)
#'
#' @export
marg_alt_g <- function(x, beta = rep(1, length(x)), lg = TRUE, ...) {
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

#' Marginal likelihood in F1 population of tetraploids
#'
#' Use this if the offspring genotypes are known, but the parental
#' genotypes are **not** known.
#'
#' @inheritParams marg_alt_g
#' @param p1_gl A vector of parent 1's genotype log-likelihoods.
#' @param p2_gl A vector of parent 2's genotype log-likelihoods.
#' @param mixprop Mixing proportion with uniform to avoid
#'     numerical issues in stan.
#'
#' @author David Gerard
#'
#' @examples
#' set.seed(1)
#' q <- hwep::zygdist(alpha = 0.1, G1 = 2, G2 = 2, ploidy = 4)
#' q <- runif(5)
#' q <- q / sum(q)
#' x <- c(stats::rmultinom(n = 1, size = 1000, prob = q))
#' p1_gl <- rep(log(0.2), 5)
#' p2_gl <- rep(log(0.2), 5)
#' mixprop <- 0.001
#' malt <- hwep:::ddirmult(x = x, alpha = rep(1, 5), lg = TRUE)
#' mnull <- marg_f1_g4(x,
#'                     chains = 1,
#'                     p1_gl = c(-10, -10, -1, -10, -10),
#'                     p2_gl = c(-10, -10, -1, -10, -10))
#' mnull - malt
#'
#' @export
marg_f1_g4 <- function(x,
                       p1_gl = rep(log(0.2), 5),
                       p2_gl = rep(log(0.2), 5),
                       mixprop = 0.001,
                       lg = TRUE, ...) {
  stopifnot(length(x) == 5,
            length(p1_gl) == 5,
            length(p2_gl) == 5,
            length(mixprop) == 1)
  stopifnot(mixprop > 0, mixprop <= 1)
  drbound <- hwep::drbounds(ploidy = 4)
  stan_dat <- list(x = x,
                   drbound = drbound,
                   p1_gl = p1_gl,
                   p2_gl = p2_gl,
                   mixprop = mixprop)
  stan_out <- rstan::sampling(object = stanmodels$marg_dr_npp_gl4,
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
#' @inheritParams marg_f1_g4
#' @param g1 The first parent's genotype.
#' @param g2 The second parent's genoype.
#'
#' @author David Gerard
#'
#' @examples
#' q <- hwep::zygdist(alpha = 0.1, G1 = 2, G2 = 2, ploidy = 4)
#' q <- runif(5)
#' q <- q / sum(q)
#' x <- c(stats::rmultinom(n = 1, size = 1000, prob = q))
#' g1 <- 2
#' g2 <- 2
#' marg_null <- marg_f1_g4_pknown(x = x, g1 = g1, g2 = g2, chains = 1)
#' marg_alt <- hwep:::ddirmult(x = x, alpha = rep(1, 5), lg = TRUE)
#' marg_null - marg_alt
#'
#' @export
marg_f1_g4_pknown <- function(x,
                              g1,
                              g2,
                              lg = TRUE, ...) {
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

#' Marginal likelihood in tetraploid F1 population
#'
#' Here, parental genotypes are known.
#'
#' @inheritParams marg_f1_g4
#' @param g1 The first parent's genotype.
#' @param g2 The second parent's genoype.
#'
#' @author David Gerard
#'
#' @examples
#' q <- hwep::zygdist(alpha = 0.1, G1 = 2, G2 = 2, ploidy = 4)
#' q <- runif(5)
#' q <- q / sum(q)
#' x <- c(stats::rmultinom(n = 1, size = 1000, prob = q))
#' g1 <- 2
#' g2 <- 2
#' marg_null <- pp_marg_f1_g4_pknown(x = x, g1 = g1, g2 = g2, chains = 1)
#' marg_alt <- hwep:::ddirmult(x = x, alpha = rep(1, 5), lg = TRUE)
#' marg_null - marg_alt
#'
#' @export
pp_marg_f1_g4_pknown <- function(x,
                              g1,
                              g2,
                              lg = TRUE, ...) {
  stopifnot(length(x) == 5,
            length(g1) == 1,
            length(g2) == 1)
  stopifnot(g1 >= 0, g1 <= 4, g2 >= 0, g2 <= 4)
  drbound <- hwep::drbounds(ploidy = 4)
  ppbound <- 1/3
  stan_dat <- list(x = x,
                   drbound = drbound,
                   ppbound = ppbound,
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
