#' Genotype frequencies of gametes when one parent's genotype is known
#'
#' @param a The double reduction rate
#' @param x The preferential pairing parameter
#' @param p The
#' @param ... Additional arguments to pass to \code{\link[rstan]{sampling}()}.
#'
#' @return The gamete genotype frequency
#'
#' @author Mira Thakkar
#'
#' @examples
#' a <- 1/6
#' x <- 1/3
#' p <- 2
#' pvec(a = a, x = x, p = p)
#'
#' @export
pvec <- function(a, x, p){
  if (p > 4 | p < 0 | is.na(p)){
    stop("Invalid input")
  } else if(p == 0){
    prob_Y0 = 1
    prob_Y1 = 0
    prob_Y2 = 0
  } else if(p == 1){
    prob_Y0 = 0.5 + (0.25 * a)
    prob_Y1 = 0.5 - (0.5 * a)
    prob_Y2 = a / 4.0
  } else if(p == 2){
    prob_Y0 = (0.5 * a) + (0.25 * (1 - a) * (1 - x))
    prob_Y1 = 0.5 * (1 - a) * (1 + x)
    prob_Y2 = (0.5 * a) + (0.25 * (1 - a) * (1 - x))
  } else if(p == 3){
    prob_Y0 = a / 4.0
    prob_Y1 = 0.5 - (0.5 * a)
    prob_Y2 = 0.5 + (0.25 * a)
  } else if(p == 4){
    prob_Y0 = 0
    prob_Y1 = 0
    prob_Y2 = 1
  }
  pv <- c(prob_Y0, prob_Y1, prob_Y2)
  return(pv)
}

#' Function that takes as input the double reduction rate, the preferential
#' pairing rate, and parent genotypes to return zygote genotype frequencies.
#'
#' @param alpha The double reduction rate
#' @param xi The preferential pairing parameter
#' @param p1 The first parent's genotype
#' @param p2 The second parent's genotype
#' @param ... Additional arguments to pass to \code{\link[rstan]{sampling}()}.
#'
#' @return Zygote genotype frequencies
#'
#' @author Mira Thakkar
#'
#' @examples
#' alpha <- 1/6
#' xi <- 1/3
#' p1 <- 2
#' p2 <- 3
#' offspring_gf(alpha = alpha, xi = xi, p1 = p1, p2 = p2)
#'
#' @export
offspring_gf <- function(alpha, xi, p1, p2){

  pvec1 <- pvec(a = alpha, x = xi, p = p1)
  pvec2 <- pvec(a = alpha, x = xi, p = p2)

  qvec <- stats::convolve(pvec1, rev(pvec2), type = "open")

  stopifnot(qvec > -1e-06)
  stopifnot(abs(sum(qvec) - 1) < 1e-06)
  qvec[qvec < 0] <- 0
  qvec <- qvec / sum(qvec)

  return(qvec)
}


#' A function called offspring_geno() Which takes as input the offspring
#' genotype frequencies and a sample size and returns simulated genotypes.
#'
#' @param x Vector of offspring genotype frequencies
#' @param n Sample size
#'
#' @return Simulated genotypes
#'
#' @author Mira Thakkar
#'
#' @examples
#' x <- offspring_gf(alpha = 1/6, xi = 1/3, p1 = 2, p2 = 3)
#' n <- 1000
#' pvec(x = x, n = n)
#'
#' @export
offspring_geno <- function(x, n){

  sim_gen <- c(stats::rmultinom(n = 1, size = n, prob = x))

  return(sim_gen)
}

#' Converts genotype counts to genotype vectors.
#'
#' @param gcount The vector of genotype counts.
#'
#' @author David Gerard
#'
#' @export
#'
#' @examples
#' gcount <- c(1, 2, 3, 0, 5)
#' gcount_to_gvec(gcount = gcount)
gcount_to_gvec <- function(gcount) {
  unlist(mapply(FUN = rep, x = seq_along(gcount) - 1, each = gcount))
}

#' Function which takes as input (i) the parent genotypes,
#' (ii) the offspring genotypes, (iii) sequencing error rate, (iv) read
#' depth, (v) bias, (vi) overdispersion and returns genotype likelihoods.
#'
#' @param genovec Offspring genotypes
#' @param p1_geno Parent 1 genotype
#' @param p2_geno Parent 2 genotype
#' @param ploidy Ploidy
#' @param seq Sequencing error rate
#' @param rd Read depth
#' @param bias Bias
#' @param od Overdispersion
#'
#' @return Genotype likelihoods
#'
#' @author Mira Thakkar
#'
#' @export
po_gl <- function(genovec, p1_geno, p2_geno, ploidy, seq = 0.01, rd = 10, bias = 1, od = 0.01) {
  n <- length(genovec)
  sizevec <- rep(rd, length.out = n)
  refvec <- updog::rflexdog(sizevec = sizevec, geno = genovec, ploidy = ploidy, seq = seq, bias = bias, od = od)
  p1ref <- updog::rflexdog(sizevec = rd, geno = p1_geno, ploidy = ploidy, seq = seq, bias = bias, od = od)
  p2ref <- updog::rflexdog(sizevec = rd, geno = p2_geno, ploidy = ploidy, seq = seq, bias = bias, od = od)

  fout <- updog::flexdog_full(refvec = refvec,
                              sizevec = sizevec,
                              ploidy = ploidy,
                              model = "f1pp",
                              seq = seq,
                              bias = bias,
                              od = od,
                              update_bias = FALSE,
                              update_seq = FALSE,
                              update_od = FALSE,
                              p1ref = p1ref,
                              p1size = rd,
                              p2ref = p2ref,
                              p2size = rd)

  return(fout)
}
