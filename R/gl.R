#############################
## Methods using genotype log-likelihoods
#############################

#' Marginal likelihood under alternative using genotype log-likelihoods
#'
#' @inheritParams marg_f1_dr_pp_gl4
#' @param beta The vector of hyperparameters.
#'
#' @return The marginal likelihood.
#'
#' @author Mira Thakkar and David Gerard
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

#' Marginal likelihood, no double reduction, no preferential pairing, genotype likelihoods.
#'
#' @inheritParams marg_f1_dr_pp_gl4
#'
#' @author Mira Thakkar and David Gerard
#'
#' @examples
#' gl <- matrix(runif(20), ncol = 5)
#' marg_f1_dr_npp_gl4(gl = gl)
#'
#' @export
marg_f1_ndr_npp_gl4 <- function(gl,
                                p1_gl = rep(log(0.2), 5),
                                p2_gl = rep(log(0.2), 5),
                                lg = TRUE) {

}

#' Marginal likelihood, double reduction, no preferential pairing, genotype likelihoods.
#'
#' @inheritParams marg_f1_dr_pp_gl4
#'
#' @author Mira Thakkar and David Gerard
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
                               lg = TRUE,
                               ...) {
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

#' Marginal likelihood, no double reduction, preferential pairing, genotype likelihoods.
#'
#' @inheritParams marg_f1_dr_pp_gl4
#'
#' @author Mira Thakkar and David Gerard
#'
#' @examples
#' gl <- matrix(runif(20), ncol = 5)
#' marg_f1_dr_npp_gl4(gl = gl)
#'
#' @export
marg_f1_ndr_pp_gl4 <- function(gl,
                               p1_gl = rep(log(0.2), 5),
                               p2_gl = rep(log(0.2), 5),
                               mixprop = 0.001,
                               lg = TRUE,
                               ...) {
  stopifnot(ncol(gl) == 5,
            length(p1_gl) == 5,
            length(p2_gl) == 5,
            length(mixprop) == 1)
  stopifnot(mixprop > 0, mixprop <= 1)
  stan_dat <- list(gl = gl,
                   N = nrow(gl),
                   p1_gl = p1_gl,
                   p2_gl = p2_gl,
                   mixprop = mixprop)
  stan_out <- rstan::sampling(object = stanmodels$marg_ndr_pp_gl4,
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

#' Marginal likelihood, double reduction, preferential pairing, genotype likelihoods.
#'
#' @param gl A matrix of genotype log-likelihoods. The rows index the
#'     individuals and the columns index the possible genotypes. So
#'     \code{gl[i, k]} is the genotype log-likelihood for individual i and
#'     genotype k-1.
#' @param p1_gl The vector of genotype likelihoods of parent 1.
#' @param p1_gl The vector of genotype likelihoods of parent 2.
#' @param mixprop The mixing proportion with the uniform for mixing purposes.
#' @param lg A logical. Should we log the marginal likelihood (\code{TRUE})
#'     or not (\code{FALSE})?
#' @param ... Additional paramters sent to \code{\link[stan]{sampling}()}.
#'
#' @author Mira Thakkar and David Gerard
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
                              lg = TRUE,
                              ...) {
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
