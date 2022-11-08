// Tetraploid F1 test

functions {
  // alpha Double reduction rate
  // g parent genotype
  // khalf ploidy / 2 + 1
  // return: gamete frequencies of a parent
  real[] segfreq4(real alpha, int g) {
    real p[3];
    if (g == 0) {
      p[1] = 1.0;
      p[2] = 0.0;
      p[3] = 0.0;
    } else if (g == 1) {
      p[1] = (2.0 + alpha) / 4.0;
      p[2] = 2.0 * (1.0 - alpha) / 4.0;
      p[3] = alpha / 4.0;
    } else if (g == 2) {
      p[1] = (1.0 + 2.0 * alpha) / 6.0;
      p[2] = 4.0 * (1.0 - alpha) / 6.0;
      p[3] = (1.0 + 2.0 * alpha) / 6.0;
    } else if (g == 3) {
      p[1] = alpha / 4.0;
      p[2] = 2.0 * (1.0 - alpha) / 4.0;
      p[3] = (2.0 + alpha) / 4.0;
    } else if (g == 4) {
      p[1] = 0.0;
      p[2] = 0.0;
      p[3] = 1.0;
    } else {
      // do nothing
    }
    return p;
  }

  // p1 gamete frequencies from parent 1.
  // p2 gamete frequencies from parent 2.
  // K ploidy
  // khalf K / 2 + 1 so stan does not complain about integer division
  vector convolve(real[] p1, real[] p2, int K, int khalf) {
    vector[K+1] q;
    for (k in 1:(K+1)) {
      int iup = min(k - 1, khalf - 1);
      int ilo = max(0, k - khalf);
      q[k] = 0.0;
      for (i in ilo:iup) {
        q[k] += p1[i + 1] * p2[k - i];
      }
    }
    return q;
  }
}

data {
  real p1_gl[5]; // genotype log-likelihoods for parent 1
  real p2_gl[5]; // genotype log-likelihoods for parent 2
  int x[5]; // genotype counts
  real drbound; // double reduction bound
}

parameters {
  real<lower=0.0,upper=drbound> alpha; // double reduction rate
}

transformed parameters {
  matrix[5, 5] glmat;
  for (i in 1:5) {
    real p1[3] = segfreq4(alpha, i);
    for (j in 1:5) {
      real p2[3] = segfreq4(alpha, i);
      vector[5] q = convolve(p1, p2, 4, 3);
      glmat[i, j] = multinomial_lpmf(x | q) + p1_gl[i] + p2_gl[j];
    }
  }

}

model {
  target += uniform_lpdf(alpha | 0.0, drbound);
  target += log_sum_exp(glmat);
}
