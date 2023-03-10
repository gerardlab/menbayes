###############
## Marginal likelihoods when parental genotypes are known, but offspring
## only have genotype likelihoods
###############

#' Marginal likelihood, no double reduction, no preferential pairing, parent genotypes known, offspring genotypes unknown
#'
#' @inheritParams marg_f1_dr_pp_glpknown4
#'
#' @examples
#' ## null sims
#' set.seed(1)
#' gf <- offspring_gf(alpha = 0, xi = 1/3, p1 = 2, p2 = 2)
#' gcount <- offspring_geno(x = gf, n = 20)
#' gvec <- gcount_to_gvec(gcount)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_ndr_npp_glpknown4(gl = gl, p1 = 2, p2 = 2)
#' malt <- -48.81 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#'
#' ## alt sims
#' set.seed(1)
#' gvec <- rep(0:4, each = 4)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_ndr_npp_glpknown4(gl = gl, p1 = 2, p2 = 2)
#' malt <- -46.82 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#'
#' @author Mira Thakkar and David Gerard
#'
#' @export
marg_f1_ndr_npp_glpknown4 <- function(gl,
                                      p1,
                                      p2,
                                      lg = TRUE) {
  gf <- offspring_gf(alpha = 0, xi = 1/3, p1 = p1, p2 = p2)
  lgf <- log(gf)
  mx <- sum(apply(X = t(gl) + lgf, MARGIN = 2, FUN = updog::log_sum_exp))
  if(!lg) {
    mx <- exp(mx)
  }
  return(mx)
}

#' Marginal likelihood, double reduction, no preferential pairing, parent genotypes known, offspring genotypes unknown
#'
#' @inheritParams marg_f1_dr_pp_glpknown4
#'
#' @examples
#' ## null sims
#' set.seed(1)
#' gf <- offspring_gf(alpha = 0, xi = 1/3, p1 = 2, p2 = 2)
#' gcount <- offspring_geno(x = gf, n = 20)
#' gvec <- gcount_to_gvec(gcount)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_dr_npp_glpknown4(gl = gl, p1 = 2, p2 = 2, chains = 1)
#' malt <- -48.81 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#'
#' ## alt sims
#' set.seed(1)
#' gvec <- rep(0:4, each = 4)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_dr_npp_glpknown4(gl = gl, p1 = 2, p2 = 2, chains = 1)
#' malt <- -46.82 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#'
#' @author Mira Thakkar and David Gerard
#'
#' @export
marg_f1_dr_npp_glpknown4 <- function(gl,
                                     p1,
                                     p2,
                                     mixprop = 0.001,
                                     lg = TRUE,
                                     ...) {
  stopifnot(ncol(gl) == 5,
            length(p1) == 1,
            length(p2) == 1)
  drbound <- hwep::drbounds(ploidy = 4)
  stan_dat <- list(gl = gl,
                   N = nrow(gl),
                   drbound = drbound,
                   g1 = p1,
                   g2 = p2,
                   mixprop = mixprop)
  stan_out <- rstan::sampling(object = stanmodels$marg_dr_npp_glpknown4,
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


#' Marginal likelihood, no double reduction, preferential pairing, parent genotypes known, offspring genotypes unknown
#'
#' @inheritParams marg_f1_dr_pp_glpknown4
#'
#' @examples
#' ## null sims
#' set.seed(1)
#' gf <- offspring_gf(alpha = 0, xi = 1/3, p1 = 2, p2 = 2)
#' gcount <- offspring_geno(x = gf, n = 20)
#' gvec <- gcount_to_gvec(gcount)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_ndr_pp_glpknown4(gl = gl, p1 = 2, p2 = 2, chains = 1)
#' malt <- -48.81 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#'
#' ## alt sims
#' set.seed(1)
#' gvec <- rep(0:4, each = 4)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_ndr_pp_glpknown4(gl = gl, p1 = 2, p2 = 2, chains = 1)
#' malt <- -46.82 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#'
#' @author Mira Thakkar and David Gerard
#'
#' @export
marg_f1_ndr_pp_glpknown4 <- function(gl,
                                     p1,
                                     p2,
                                     mixprop = 0.001,
                                     lg = TRUE,
                                     ...) {
  stopifnot(ncol(gl) == 5,
            length(p1) == 1,
            length(p2) == 1)
  stan_dat <- list(gl = gl,
                   N = nrow(gl),
                   g1 = p1,
                   g2 = p2,
                   mixprop = mixprop)
  stan_out <- rstan::sampling(object = stanmodels$marg_ndr_pp_glpknown4,
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

#' Marginal likelihood, double reduction, preferential pairing, parent genotypes known, offspring genotypes unknown
#'
#' @param gl A matrix of offspring genotype log-likelihoods. The rows index the
#'     individuals and the columns index the possible genotypes. So
#'     \code{gl[i, k]} is the offspring genotype log-likelihood for individual i and
#'     genotype k-1.
#' @param p1 Genotype of parent 1.
#' @param p2 Genotype of parent 2.
#' @param mixprop The mixing proportion with the uniform for mixing purposes.
#' @param lg A logical. Should we log the marginal likelihood (\code{TRUE}) or
#'     not (\code{FALSE})?
#' @param ... Additional paramters sent to \code{\link[stan]{sampling}()}.
#'
#' @examples
#' ## null sims
#' set.seed(1)
#' gf <- offspring_gf(alpha = 0, xi = 1/3, p1 = 2, p2 = 2)
#' gcount <- offspring_geno(x = gf, n = 20)
#' gvec <- gcount_to_gvec(gcount)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_dr_pp_glpknown4(gl = gl, p1 = 2, p2 = 2, chains = 1)
#' malt <- -48.81 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#'
#' ## alt sims
#' set.seed(1)
#' gvec <- rep(0:4, each = 4)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_dr_pp_glpknown4(gl = gl, p1 = 2, p2 = 2, chains = 1)
#' malt <- -46.82 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#'
#' @author Mira Thakkar and David Gerard
#'
#' @export
marg_f1_dr_pp_glpknown4 <- function(gl,
                                    p1,
                                    p2,
                                    mixprop = 0.001,
                                    lg = TRUE,
                                    ...) {
  stopifnot(ncol(gl) == 5,
            length(p1) == 1,
            length(p2) == 1)
  drbound <- hwep::drbounds(ploidy = 4)
  stan_dat <- list(gl = gl,
                   N = nrow(gl),
                   drbound = drbound,
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
