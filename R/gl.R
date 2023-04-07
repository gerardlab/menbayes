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
#' gl <- hwep::simgl(nvec = x, rdepth = 100, od = 0, bias = 1, seq = 0)
#' marg_alt_gl(gl = gl, chains = 1)
#'
#' @export
marg_alt_gl <- function(gl, beta = rep(1, 5), lg = TRUE, ...) {
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
#' ## null sims
#' set.seed(1)
#' gf <- offspring_gf(alpha = 0, xi = 1/3, p1 = 2, p2 = 2)
#' gcount <- offspring_geno(x = gf, n = 20)
#' gvec <- gcount_to_gvec(gcount)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_ndr_npp_gl4(gl = gl, p1_gl = c(-Inf, -Inf, 0, -Inf, -Inf), p2_gl = c(-Inf, -Inf, 0, -Inf, -Inf))
#' malt <- -48.81 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#'
#' ## alt sims
#' set.seed(1)
#' gvec <- rep(0:4, each = 4)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_ndr_npp_gl4(gl = gl, p1_gl = c(-Inf, -Inf, 0, -Inf, -Inf), p2_gl = c(-Inf, -Inf, 0, -Inf, -Inf))
#' malt <- -46.82 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#'
#' @export
marg_f1_ndr_npp_gl4 <- function(gl,
                                p1_gl = rep(log(0.2), 5),
                                p2_gl = rep(log(0.2), 5),
                                lg = TRUE) {
  stopifnot(ncol(gl) == 5,
            length(p1_gl) == 5,
            length(p2_gl) == 5)
  p1_gl <- p1_gl - updog::log_sum_exp(p1_gl)
  p2_gl <- p2_gl - updog::log_sum_exp(p2_gl)
  glmat <- matrix(data = NA_real_, nrow = 5, ncol = 5)
  for (i in 0:4) {
    for (j in 0:4) {
      gf <- offspring_gf(alpha = 0, xi = 1/3, p1 = i, p2 = j)
      lgf <- log(gf)
      glmat[i + 1, j + 1] <- sum(apply(X = t(gl) + lgf, MARGIN = 2, FUN = updog::log_sum_exp)) +
        p1_gl[[i + 1]] +
        p2_gl[[j + 1]]
    }
  }
  mx <- updog::log_sum_exp(glmat)
  if(!lg) {
    mx <- exp(mx)
  }
  return(mx)
}

#' Marginal likelihood, double reduction, no preferential pairing, genotype likelihoods.
#'
#' @inheritParams marg_f1_dr_pp_gl4
#'
#' @author Mira Thakkar and David Gerard
#'
#' @examples
#' \dontrun{
#' ## null sims
#' set.seed(1)
#' gf <- offspring_gf(alpha = 0, xi = 1/3, p1 = 2, p2 = 2)
#' gcount <- offspring_geno(x = gf, n = 20)
#' gvec <- gcount_to_gvec(gcount)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_dr_npp_gl4(gl = gl, p1_gl = c(-Inf, -Inf, 0, -Inf, -Inf), p2_gl = c(-Inf, -Inf, 0, -Inf, -Inf), chains = 1)
#' malt <- -48.81 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#'
#' ## alt sims
#' set.seed(1)
#' gvec <- rep(0:4, each = 4)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_dr_npp_gl4(gl = gl, p1_gl = c(-Inf, -Inf, 0, -Inf, -Inf), p2_gl = c(-Inf, -Inf, 0, -Inf, -Inf), chains = 1)
#' malt <- -46.82 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#' }
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
  p1_gl <- p1_gl - updog::log_sum_exp(p1_gl)
  p2_gl <- p2_gl - updog::log_sum_exp(p2_gl)
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
#' \dontrun{
#' ## null sims
#' set.seed(1)
#' gf <- offspring_gf(alpha = 0, xi = 1/3, p1 = 2, p2 = 2)
#' gcount <- offspring_geno(x = gf, n = 20)
#' gvec <- gcount_to_gvec(gcount)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_ndr_pp_gl4(gl = gl, p1_gl = c(-Inf, -Inf, 0, -Inf, -Inf), p2_gl = c(-Inf, -Inf, 0, -Inf, -Inf), chains = 1)
#' malt <- -48.81 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#'
#' ## alt sims
#' set.seed(1)
#' gvec <- rep(0:4, each = 4)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_ndr_pp_gl4(gl = gl, p1_gl = c(-Inf, -Inf, 0, -Inf, -Inf), p2_gl = c(-Inf, -Inf, 0, -Inf, -Inf), chains = 1)
#' malt <- -46.82 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#' }
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
  p1_gl <- p1_gl - updog::log_sum_exp(p1_gl)
  p2_gl <- p2_gl - updog::log_sum_exp(p2_gl)
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
#' @param ... Additional parameters sent to \code{\link[stan]{sampling}()}.
#'
#' @author Mira Thakkar and David Gerard
#'
#' @examples
#' \dontrun{
#' ## null sims
#' set.seed(1)
#' gf <- offspring_gf(alpha = 0, xi = 1/3, p1 = 2, p2 = 2)
#' gcount <- offspring_geno(x = gf, n = 20)
#' gvec <- gcount_to_gvec(gcount)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_dr_pp_gl4(gl = gl, p1_gl = c(-Inf, -Inf, 0, -Inf, -Inf), p2_gl = c(-Inf, -Inf, 0, -Inf, -Inf), chains = 1)
#' malt <- -48.81 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#'
#' ## alt sims
#' set.seed(1)
#' gvec <- rep(0:4, each = 4)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_dr_pp_gl4(gl = gl, p1_gl = c(-Inf, -Inf, 0, -Inf, -Inf), p2_gl = c(-Inf, -Inf, 0, -Inf, -Inf), chains = 1)
#' malt <- -46.82 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#' }
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
  p1_gl <- p1_gl - updog::log_sum_exp(p1_gl)
  p2_gl <- p2_gl - updog::log_sum_exp(p2_gl)
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
