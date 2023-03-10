#############################
## Methods using genotype log-likelihoods
#############################

#' Marginal likelihood under alternative using genotype log-likelihoods
#'
#' @inheritParams marg_alt_gl
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

#' Marginal likelihood under null using genotype log-likelihoods
#'
#' @inheritParams marg_f1_dr_npp_gl4
#' @param gl A matrix of genotype log-likelihoods. The rows index the
#'     individuals and the columns index the possible genotypes. So
#'     \code{gl[i, k]} is the genotype log-likelihood for individual i and
#'     genotype k-1.
#'
#' @examples
#' gl <- matrix(runif(20), ncol = 5)
#' marg_f1_dr_npp_gl4(gl = gl)
#'
#' @export
marg_f1_dr_npp_gl4 <- function(gl,
                        p1_gl = rep(log(0.2), 5),
                        p2_gl = rep(log(0.2), 5),
                        mixprop = 0.001,
                        lg = TRUE, ...) {
  stopifnot(ncol(gl) == 5,
            length(p1_gl) == 5,
            length(p2_gl) == 5,
            length(mixprop) == 1)
  stopifnot(mixprop > 0, mixprop <= 1)
  drbound <- hwep::drbounds(ploidy = 4)
  stan_dat <- list(gl = gl,
                   N = nrow(gl),
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

#' Marginal likelihood under null using genotype log-likelihoods
#'
#' @inheritParams marg_f1_dr_npp_gl4
#' @param gl A matrix of genotype log-likelihoods. The rows index the
#'     individuals and the columns index the possible genotypes. So
#'     \code{gl[i, k]} is the genotype log-likelihood for individual i and
#'     genotype k-1.
#'
#' @examples
#' gl <- matrix(runif(20), ncol = 5)
#' marg_f1_dr_pp_gl4(gl = gl)
#'
#' @export
marg_f1_dr_pp_gl4 <- function(gl,
                        p1_gl = rep(log(0.2), 5),
                        p2_gl = rep(log(0.2), 5),
                        mixprop = 0.001,
                        lg = TRUE, ...) {
  stopifnot(ncol(gl) == 5,
            length(p1_gl) == 5,
            length(p2_gl) == 5,
            length(mixprop) == 1)
  stopifnot(mixprop > 0, mixprop <= 1)
  drbound <- hwep::drbounds(ploidy = 4)
  ppbound <- 1/3
  stan_dat <- list(gl = gl,
                   N = nrow(gl),
                   drbound = drbound,
                   ppbound = ppbound,
                   p1_gl = p1_gl,
                   p2_gl = p2_gl,
                   mixprop = mixprop)
  stan_out <- rstan::sampling(object = stanmodels$marg_dr_pp_gl4,
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
                              mixprop = 0.001,
                              lg = TRUE, ...) {
  stopifnot(ncol(gl) == 5,
            length(p1) == 1,
            length(p2) == 1,
            length(mixprop) == 1)
  stopifnot(mixprop > 0, mixprop <= 1)
  drbound <- hwep::drbounds(ploidy = 4)
  ppbound <- 1/3
  stan_dat <- list(gl = gl,
                   N = nrow(gl),
                   drbound = drbound,
                   ppbound = ppbound,
                   g1 = p1,
                   g2 = p2,
                   mixprop = mixprop)
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
