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
#' @return The mariginal likelihood.
#'
#' @author David Gerard
#'
#' @examples
#' x <- c(1, 2, 7, 3, 3)
#' beta <- rep(1, 5)
#' marg_alt_g(x = x, beta = beta, chains = 1)
#'
#' @export
marg_alt_g <- function(x, beta, lg = TRUE, ...) {
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
#' @inheritParams marg_alt_g
#' @param p1_gl A vector of parent 1's genotype log-likelihoods.
#' @param p2_gl A vector of parent 2's genotype log-likelihoods
#'
#' @author David Gerard
#'
#' @examples
#' set.seed(1)
#' q <- hwep::zygdist(alpha = 0.1, G1 = 2, G2 = 2, ploidy = 4)
#' x <- c(stats::rmultinom(n = 1, size = 10, prob = q))
#' p1_gl <- rep(log(0.2), 5)
#' p2_gl <- rep(log(0.2), 5)
#'
#' @export
marg_f1_g4 <- function(x,
                       p1_gl = rep(log(0.2), 5),
                       p2_gl = rep(log(0.2), 5),
                       lg = TRUE, ...) {
  stopifnot(length(x) == 5,
            length(p1_gl) == 5,
            length(p2_gl) == 5)
  drbound <- hwep::drbounds(ploidy = 4)
  stan_dat <- list(x = x, drbound = drbound, p1_gl = p1_gl, p2_gl = p2_gl)
  stan_out <- rstan::sampling(object = stanmodels$f1_g4,
                              data = stan_dat,
                              verbose = FALSE,
                              show_messages = FALSE,
                              init = 0)
  bridge_out <- bridgesampling::bridge_sampler(stan_out, verbose = FALSE, silent = TRUE)

  if (lg) {
    mx <- bridge_out$logml
  } else {
    mx <- exp(bridge_out$logml)
  }
  return(mx)
}
