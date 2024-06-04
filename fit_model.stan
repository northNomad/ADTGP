data{
  int N;
  array[N] real<lower=0> igg;
  array[N] int poi;

  int nrep;
  array[nrep] real<lower=0> igg_rep;
}

parameters{
  real mu0;
  real<lower=0> thi;

  real<lower=0> etasq;
  real<lower=0> rhosq;
  real<lower=0> sigmasq;
}

transformed parameters{
  matrix[N, N] K;
  K = gp_exp_quad_cov(igg, sqrt(etasq), rhosq);
  for(i in 1:N){
    K[i, i] = K[i, i] + sigmasq;
  }
}

model{
  mu0 ~ gamma(7, 2);
  thi ~ gamma(.5, .5);
  etasq ~ normal(2, 1);
  rhosq ~ normal(2, 1);
  sigmasq ~ normal(2, 1);

  vector[N] gamma;
  gamma ~ multi_normal_cholesky(rep_vector(0, N), cholesky_decompose(K));

  for(i in 1:N){
    poi[i] ~ neg_binomial_2(exp(mu0 + gamma[i]), thi);
  }
}

generated quantities{
  array[nrep] int poi_given_igg;
  matrix[nrep, nrep] K_given_igg;
  K_given_igg = gp_exp_quad_cov(igg_rep, sqrt(etasq), rhosq);
  for(i in 1:nrep){
    K_given_igg[i ,i] = K_given_igg[i, i] + sigmasq;
  }
  vector[nrep] gamma_rep;
  gamma_rep = multi_normal_cholesky_rng(rep_vector(0, nrep),
                                        cholesky_decompose(K_given_igg));
  for(i in 1:nrep){
    poi_given_igg[i] = neg_binomial_2_rng(exp(mu0 + gamma_rep[i]), thi);
  }
}
