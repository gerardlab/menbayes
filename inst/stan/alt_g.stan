data {
  int<lower=2> K; // ploidy
  int x[K+1]; // genotype counts
  vector<lower=0>[K+1] beta; // prior concentration parameters
}

parameters {
  simplex[K + 1] q;
}

model {
  target += dirichlet_lpdf(q | beta); // **CANNOT** use "~" notation because drops terms, need to use target +=
  target += multinomial_lpmf(x | q);
}
