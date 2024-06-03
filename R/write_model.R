ADTnorm_prior_predict <- function(
    n=NULL,
    igg=NULL,
    priors=list(mu0="gamma(7, 2)",
                thi="normal(.5, .5)",
                etasq="normal(2, 2)",
                rhosq="normal(2, 2)",
                sigmasq="normal(0, 1)"
           ),
    fn="StanModel_prior_predict.stan",
    parallel_chains=4,
    iter_warmup=1e3,
    iter_sampling=1e3,
    variables="poi"
    ){

      ## Checks if 'n' is type integer
      if(is.integer(n) == FALSE){ stop("Error: n is not an integer.") }

      ## Checks if 'igg' is type real
      if(is.numeric(igg) == FALSE){ stop("Error: igg is not numeric.") }

      ## Checks if length of 'n' and 'igg' matches
      if(n != length(igg)){ stop("Error: length of igg and n doesn't match.")}

      priors <- lapply(priors, function(i) paste0("~", i, ";")) ## change syntax
      priors <- lapply(seq_along(priors),
                       function(y, n, i){ paste(n[[i]], y[[i]]) },
                       y=priors, n=names(priors)
                       )
      priors <- unlist(priors)

      c("data{
          int N;
          array[N] real igg;
        }",

        "",

        "parameters{
          real mu0;
          real<lower=0> thi;

          real<lower=0> etasq; //GP priors
          real<lower=0> rhosq;
          real<lower=0> sigmasq;
        }",

        "",

        "transformed parameters{
          matrix[N, N] K;
          K = gp_exp_quad_cov(igg, sqrt(etasq), rhosq);
          for(i in 1:N){
            K[i, i] = K[i, i] + sigmasq;
          }
        }",

        "",

        "model{", priors, "}",

        "",

        "generated quantities{
          vector[N] poi;
          vector[N] gamma;
          gamma = multi_normal_rng(rep_vector(0, N), K);

          for(i in 1:N){
            poi[i] = neg_binomial_2_log_rng(mu0 + gamma[i], thi);
          }
        }",
        ""
        ) -> model_text

      writeLines(model_text, fn) ## write model to stan file

      m <- cmdstanr::cmdstan_model(fn) ## compile model
      message("Stan model finished compiling.")

      p <- m$sample(data=list(N=n, igg=igg),
                              parallel_chains=parallel_chains,
                              iter_warmup=iter_warmup,
                              iter_sampling=iter_sampling
                    )

      list(stanfit=p,
           draws=p$draws(variables=variables, format="draws_matrix")
           ) -> out

      return(out)
}
