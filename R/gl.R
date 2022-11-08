#############################
## Methods using genotype log-likelihoods
#############################

#' Marginal likelihood under alternative using gentoype log-likelihoods
#'
#' @inheritParams marg_alt_g
#' @param gl A matrix of genotype log-likelihoods. The rows index the
#'     individuals and the columns index the possible genotypes. So
#'     \code{gl[i, k]} is the genotype log-likelihood for individual i and
#'     genotype k-1.
#'
#' @return The marginal likelihood.
#'
#' @author David Gerard
#'
#' @examples
#' x <- c(1, 2, 7, 3, 3)
#' beta <- rep(1, 5)
#' gl <- hwep::simgl(nvec = x, rdepth = 100, od = 0, bias = 1, seq = 0)
#' marg_alt_gl(gl = gl, beta = beta, chains = 1)
#'
#' @export
marg_alt_gl <- function(gl, beta, lg = TRUE, ...) {
  ploidy <- ncol(gl) - 1
  nind <- nrow(gl)
  stopifnot(ncol(gl) == length(beta))
  stan_dat <- list(N = nind, K = ploidy, gl = gl, beta = beta)
  stan_out <- rstan::sampling(object = stanmodels$alt_gl,
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
