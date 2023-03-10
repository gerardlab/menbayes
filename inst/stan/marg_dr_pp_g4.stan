// Tetraploid F1 test
// Double reduction and preferential pairing
// parental genotypes are known

functions {
  // alpha Double reduction rate
  // xi Preferential pairing rate
  // g parent genotype
  // khalf ploidy / 2 + 1
  // return: gamete frequencies of a parent
  vector segfreq4(real alpha, real xi, int g) {
    vector[3] p;
    if (g == 0) {
      p[1] = 1.0;
      p[2] = 0.0;
      p[3] = 0.0;
    } else if (g == 1) {
      p[1] = 0.5 + (0.25 * alpha);
      p[2] = 0.5 - (0.5 * alpha);
      p[3] = alpha / 4.0;
    } else if (g == 2) {
      p[1] = (0.5 * alpha) + (0.25 * (1 - alpha) * (1 - xi));
      p[2] = 0.5 * (1 - alpha) * (1 + xi);
      p[3] = (0.5 * alpha) + (0.25 * (1 - alpha) * (1 - xi));
    } else if (g == 3) {
      p[1] = alpha / 4.0;
      p[2] = 0.5 - (0.5 * alpha);
      p[3] = 0.5 + (0.25 * alpha);
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
  vector convolve(vector p1, vector p2, int K, int khalf) {
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
  int<lower=0> x[5]; // genotype counts
  real<lower=0.0,upper=1.0> drbound; // upper bound of double reduction rate
  int<lower=0,upper=4> g1; // first parent genotype
  int<lower=0,upper=4> g2; // second parent genotype
}

parameters {
  real<lower=0,upper=drbound> alpha; // double reduction rate
  real<lower=0,upper=1> xi; // preferential pairing rate
}

model {
  vector[3] p1;
  vector[3] p2;
  vector[5] q;
  p1 = segfreq4(alpha, xi, g1);
  p2 = segfreq4(alpha, xi, g2);
  q = convolve(p1, p2, 4, 3);
  target += uniform_lpdf(alpha | 0.0, drbound);
  target += beta_lpdf(xi | 1.0, 2.0);
  target += multinomial_lpmf(x | q);
}
