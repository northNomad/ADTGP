ADTnorm_lfe <- function(
    igg=NULL,
    poi=NULL,
    design_matrix = NULL,
    fn="StanModel_fit_data.stan",
    parallel_chains=4,
    iter_warmup=5e3,
    iter_sampling=5e3){

  N <- length(igg)

  covariates <- colnames(design_matrix)

  ## checks if first column of design_matrix is mu0
  if(colnames(design_matrix)[1] != "mu0"){
    stop("Error: first column of design_matrix is not mu0")
  }

  ## Check how many levels are in each covariate
  n_levels <- apply(design_matrix, 2, function(x) length(unique(x)))

  ## Create parameters for each covariate to pass to stan
  b_covariates <- paste0("b_", covariates[-1])
  stan_parameters <- paste0("vector[", n_levels[-1], "] ", b_covariates, ";")

  ## Create priors for each covariate to pass to stan
  stan_priors <- paste0(b_covariates, " ~ std_normal();")

  ## Create generated quantities block
  covariates_rep <- paste0(covariates[-1], "_rep")
  data_rep <- paste0("vector[n_rep] ", covariates_rep, ";")
  cov_lfe <- paste0(b_covariates, "[", covariates_rep, "[i]]", collapse = " + ")

  c("generated quantities{
    vector[n_rep] poi_rep;
    vector[n_rep] mu_rep;
    vector[n_rep] gamma_rep;
    vector[n_rep] theta_rep;
    matrix[n_rep, n_rep] K_rep;
    matrix[n_rep, n_rep] LK_rep;

    K_rep = gp_exp_quad_cov(igg_rep, sqrt(etasq), sqrt(rhosq));
    for(i in 1:n_rep){
      K_rep[i, i] = K_rep[i, i] + delta;
      theta_rep[i] = std_normal_rng();
    }
    LK_rep = cholesky_decompose(K_rep);
    gamma_rep = LK_rep * theta_rep;",

    "for(i in 1:n_rep){
      mu_rep[i] = mu0 + gamma_rep[i] +", cov_lfe, ";",
    "}",
    "for(i in 1:n_rep){
     poi_rep[i] = neg_binomial_2_log_rng(mu_rep[i], phi);
   }",
    "}") -> stan_block_generated_quantities


  ## Create data block
  stan_covariates_data_block <- paste0("array[N] int ", covariates[-1], ";")
  stan_covariates_rep_data_block <- paste0("array[n_rep] int ", covariates_rep, ";")
  c("data{
    int N;
    int n_rep;
    array[n_rep] real igg_rep;
    array[N] real igg;
    array[N] int poi;",
    stan_covariates_data_block,
    stan_covariates_rep_data_block,
    "}") -> stan_block_data

  ## Create transformed data block
  "transformed data{
    real delta = 1e-9;
  }" -> stan_block_transformed_data

  ## Create parameters block
  c("parameters{
    real mu0;
    real<lower=0> phi;
    real<lower=0> etasq; //GP priors
    real<lower=0> rhosq;
    vector[N] theta; // isotropic normal",
    stan_parameters,
  "}") -> stan_block_parameters

  ## Create transformed parameters block
  "transformed parameters{
    matrix[N, N] K;
    matrix[N, N] L_K;
    K = gp_exp_quad_cov(igg, sqrt(etasq), sqrt(rhosq));
    for(i in 1:N){
      K[i, i] = K[i, i] + delta;
    }
    L_K = cholesky_decompose(K);
  }" -> stan_block_transformed_parameters

  ## Create model block
  stan_lfe <- paste0("mu[i] = mu0 + gamma[i] + ",
                     paste0(b_covariates, "[", covariates[-1], "[i]]",
                            collapse = " + "), ";")
  c("model{
    vector[N] gamma;
    theta ~ std_normal();
    gamma = L_K * theta;

    mu0 ~ gamma(7, 2);
    phi ~ gamma(3, 1);
    etasq ~ normal(2, 1);
    rhosq ~ uniform(0, 5000);",
    stan_priors,

    "vector[N] mu;
    for(i in 1:N){",
    stan_lfe,
    "}",
    "poi ~ neg_binomial_2_log(mu, phi);",
  "}") -> stan_block_model



  model_text <- c(stan_block_data,
                  stan_block_transformed_data,
                  stan_block_parameters,
                  stan_block_transformed_parameters,
                  stan_block_model,
                  stan_block_generated_quantities)

  writeLines(model_text, fn) ## write model to stan file

  m <- cmdstanr::cmdstan_model(fn) ## compile model
  message("Stan model finished compiling.")

  ## Prepare data
  data <- list(N=N, poi=poi, igg=igg)
  data_covariates <- lapply(seq_len(ncol(design_matrix[, -1, drop=FALSE])),
                            function(i){
                              design_matrix[, -1, drop=FALSE][,i]
                            })

  n_rep <- prod(n_levels[-1])
  igg_rep <- rep(0, n_rep)
  data_generated_quantities <- list(n_rep=n_rep, igg_rep=igg_rep)
  data_generated_quantities2 <- expand.grid(apply(design_matrix, 2, unique)[-1])
  lapply(seq_len(ncol(data_generated_quantities2)),
         function(i){
           data_generated_quantities2[,i]
         }) -> data_generated_quantities2
  names(data_generated_quantities2) <- covariates_rep

  data_generated_quantities <- c(data_generated_quantities,
                                 data_generated_quantities2)

  names(data_covariates) <- covariates[-1]
  data <- c(data, data_covariates, data_generated_quantities)

  p <- m$sample(data=data,
                parallel_chains=parallel_chains,
                iter_warmup=iter_warmup,
                iter_sampling=iter_sampling
  )


  draws <- p$draws(variables="poi_rep", format="draws_matrix")

  poi_rep_name <- Map(paste,
                      names(data_generated_quantities2),
                      data_generated_quantities2, sep = "_")
  poi_rep_name <- do.call(paste,
                          c(data.frame(poi_rep_name), sep="_"))

  colnames(draws) <- poi_rep_name
  out <- list(stanfit=p, draws=draws)
  return(out)
}
