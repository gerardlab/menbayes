###################
## Methods when genotypes are known
###################

#' Marginal likelihood under alternative when genotypes are known
#'
#' @param x The genotype counts
#' @param alpha The prior hyperparameters
#' @param lg A logical. Should we log the marginal likelihood or not?
#' @param ... Additional arguments to pass to \code{\link[rstan]{sampling}()}.
#'
#' @return The mariginal likelihood.
#'
#' @author David Gerard
#' @export
marg_alt_g <- function(x, alpha, lg = TRUE, ...) {
  ploidy <- length(x) - 1
  stopifnot(length(x) == length(alpha))
  stan_dat <- list(K = ploidy, x = x, alpha = alpha)
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
