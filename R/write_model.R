#' Perform prior predictive check of protein count
#'
#' Writes a stan model and draws samples from the prior distribution of
#' protein count given a vector of isotype control (igg) count.
#'
#' @param n An integer declaring the number of cells to simulate.
#' @param igg A numeric vector of length \code{n} giving the igg count.
#' @param priors Likelihood functions in stan syntax as a named
#' list of characters. See reference for avaiable options.
#' \describe{
#'   \item{\strong{mu0: the grand mean in log scale}}
#'   {Default is "gamma(7, 2)"}
#'   \item{\strong{thi: overdispersion}}
#'   {Default is "gamma(.5, .5)"}
#'   \item{\strong{etasq: maximum covariance between any two cells}}
#'   {Default is "normal(2, 1)"}
#'   \item{\strong{rhosq: rate of decline in covariance}}
#'   {Default is "normal(2, 1)"}
#'   \item{\strong{sigmasq: extra covariance of the same cell}}
#'   {Default is "normal(2, 1)"}
#' }
#' @param fn File path to store stan model file.
#' @param parallel_chains Number of parallel chains to run. Defaults to 4.
#' @param iter_warmup Number of warmup iterations before sampling.
#'                    Defaults to 1,000.
#' @param iter_sampling Number of samples to draw from each chain.
#'                    Defaults to 1,000.
#' @param variables A character vector declaring which sampled variables
#'                  to return.
#'                  Defaults to "poi" which is the count of the
#'                  protein of interest.
#' @return A stanfit object and a matrix of draws stored in a named list.
#' @references \url{https://mc-stan.org/docs/reference-manual/}
#' @examples
#' \dontrun{
#' n <- 100L
#' igg <- rnbinom(100, mu=50, size=2)
#'
#' prior_samples <- ADTnorm_prior_predict(n=n, igg=igg, variables="poi")
#' mu_poi <- apply(prior_samples$draws, 2, mean)
#' plot(igg, mu_poi)
#' }
ADTnorm_prior_predict <- function(
    n=NULL,
    igg=NULL,
    priors=list(mu0="gamma(7, 2)",
                thi="gamma(.5, .5)",
                etasq="normal(2, 1)",
                rhosq="normal(2, 1)",
                sigmasq="normal(2, 1)"
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
