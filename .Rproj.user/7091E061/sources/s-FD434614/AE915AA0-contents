functions {
  // nothing for this example
}
data{
  int<lower=1> N;
  int<lower=1> p;
  matrix[N,p] X;
  vector[N] y;
  vector[p] m0;
  matrix[p,p] prec0;
  real<lower=0> nu0;
  real<lower=0> sig20;
}

transformed data{
  // nothing for this example
}

parameters{
  vector[p] beta;
  real<lower=0> sigma2;
}

transformed parameters{
  real<lower=0> sigma=sqrt(sigma2);
}

model{
  sigma2 ~ scaled_inv_chi_square(nu0, sig20);
  y ~ normal(X*beta, sigma);
  beta ~ multi_normal_prec(m0, prec0);
}

generated quantities{
  // nothing for this example
}
