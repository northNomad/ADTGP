#' Fits model and draws posterior samples of ADT counts
#'
#' Writes a stan model and draws samples from the posterior distribution of
#' protein count given a vector of isotype control (igg) count.
#' Used when only estimating the grand mean "mu0" (i.e., does not support design matrix).
#'
#' @param igg A numeric vector giving the igg count.
#' @param poi A numeric vector giving the protein of interest count.
#' @param priors Likelihood functions in stan syntax as a named
#' list of characters. See reference for avaiable options.
#' \describe{
#'   \item{\strong{mu0: the grand mean in log scale}}
#'   {Default is "gamma(7, 2)"}
#'   \item{\strong{phi: overdispersion}}
#'   {Default is "gamma(.5, .5)"}
#'   \item{\strong{etasq: maximum covariance between any two cells}}
#'   {Default is "normal(2, 1)"}
#'   \item{\strong{rhosq: rate of decline in covariance}}
#'   {Default is "uniform(0, 5000)"}
#' }
#' @param fn File path to store stan model file.
#' @param parallel_chains Number of parallel chains to run. Defaults to 4.
#' @param iter_warmup Number of warmup iterations before sampling.
#'                    Defaults to 5,000.
#' @param iter_sampling Number of samples to draw from each chain.
#'                    Defaults to 5,000.
#' @param variables A character vector declaring which sampled variables
#'                  to return.
#'                  Defaults to "poi_rep" which is the count of the
#'                  protein of interest.
#' @return A stanfit object and a matrix of draws stored in a named list.
#' @references \url{https://mc-stan.org/docs/reference-manual/}
ADTGP_InterceptOnly <- function(
    igg=NULL,
    poi=NULL,
    design_matrix = NULL,
    priors=list(mu0="gamma(7, 2)",
                phi="gamma(.5, .5)",
                etasq="normal(2, 1)",
                rhosq="uniform(0, 5000)"
    ),
    fn="StanModel_fit_data.stan",
    parallel_chains=4,
    iter_warmup=5e3,
    iter_sampling=5e3){

  N <- length(igg)

  ## Create priors for each covariate to pass to stan

  priors <- lapply(priors, function(i) paste0("~", i, ";")) ## change syntax
  priors <- lapply(seq_along(priors),
                   function(y, n, i){ paste(n[[i]], y[[i]]) },
                   y=priors, n=names(priors)
  )
  priors <- unlist(priors)

  ## Create generated quantities block
  # covariates_rep <- paste0(covariates[-1], "_rep")
  # data_rep <- paste0("vector[n_rep] ", covariates_rep, ";")
  # cov_lfe <- paste0(b_covariates, "[", covariates_rep, "[i]]", collapse = " + ")

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
      mu_rep[i] = mu0 + gamma_rep[i];",
    "}",
    "for(i in 1:n_rep){
     poi_rep[i] = neg_binomial_2_log_rng(mu_rep[i], phi);
   }",
    "}") -> stan_block_generated_quantities


  ## Create data block
  c("data{
    int N;
    int n_rep;
    array[n_rep] real igg_rep;
    array[N] real igg;
    array[N] int poi;",
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
  stan_lfe <- paste0("mu[i] = mu0 + gamma[i];")
  c("model{
    vector[N] gamma;
    theta ~ std_normal();
    gamma = L_K * theta;",
    priors,
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

  n_rep <- 10
  igg_rep <- rep(0, n_rep)
  data_generated_quantities <- list(n_rep=n_rep, igg_rep=igg_rep)

  data <- c(data, data_generated_quantities)

  p <- m$sample(data=data,
                parallel_chains=parallel_chains,
                iter_warmup=iter_warmup,
                iter_sampling=iter_sampling
  )


  draws <- p$draws(variables="poi_rep", format="draws_matrix")

  out <- list(stanfit=p, draws=draws)
  return(out)
}
