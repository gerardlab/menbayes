###############
## Marginal likelihoods when parental genotypes are known, but offspring
## only have genotype likelihoods
###############


marg_f1_ndr_npp_glpknown4 <- function(p1,
                                      p2,
                                      lg = TRUE,
                                      ...) {

}

marg_f1_dr_npp_glpknown4 <- function(p1,
                                     p2,
                                     lg = TRUE,
                                     ...) {

}

marg_f1_ndr_pp_glpknown4 <- function(p1,
                                     p2,
                                     lg = TRUE,
                                     ...) {

}

#' Marginal likelihood under null using offspring genotype log-likelihoods and parent genotypes
#'
#' @param gl A matrix of offspring genotype log-likelihoods. The rows index the
#'     individuals and the columns index the possible genotypes. So
#'     \code{gl[i, k]} is the offspring genotype log-likelihood for individual i and
#'     genotype k-1.
#' @param p1 Genotype of parent 1.
#' @param p2 Genotype of parent 2.
#'
#' @examples
#' gl <- matrix(runif(20), ncol = 5)
#' marg_f1_dr_pp_gl4(gl = gl, p1 = 2, p2 = 3)
#'
#' @export
marg_f1_dr_pp_glpknown4 <- function(gl,
                                    p1,
                                    p2,
                                    lg = TRUE,
                                    ...) {
  stopifnot(ncol(gl) == 5,
            length(p1) == 1,
            length(p2) == 1)
  drbound <- hwep::drbounds(ploidy = 4)
  ppbound <- 1/3
  stan_dat <- list(gl = gl,
                   N = nrow(gl),
                   drbound = drbound,
                   ppbound = ppbound,
                   g1 = p1,
                   g2 = p2)
  stan_out <- rstan::sampling(object = stanmodels$marg_dr_pp_glpknown4,
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
