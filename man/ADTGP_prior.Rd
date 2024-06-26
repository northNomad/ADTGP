% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Prior.R
\name{ADTGP_prior}
\alias{ADTGP_prior}
\title{Perform prior predictive check of protein count}
\usage{
ADTGP_prior(
  n = NULL,
  igg = NULL,
  priors = list(mu0 = "gamma(7, 2)", phi = "gamma(.5, .5)", etasq = "normal(2, 1)", rhosq
    = "uniform(0, 5000)"),
  fn = "StanModel_prior_predict.stan",
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  variables = "poi"
)
}
\arguments{
\item{n}{An integer declaring the number of cells to simulate.}

\item{igg}{A numeric vector of length \code{n} giving the igg count.}

\item{priors}{Likelihood functions in stan syntax as a named
list of characters. See reference for avaiable options.
\describe{
  \item{\strong{mu0: the grand mean in log scale}}
  {Default is "gamma(7, 2)"}
  \item{\strong{phi: overdispersion}}
  {Default is "gamma(.5, .5)"}
  \item{\strong{etasq: maximum covariance between any two cells}}
  {Default is "normal(2, 1)"}
  \item{\strong{rhosq: rate of decline in covariance}}
  {Default is "uniform(0, 5000)"}
}}

\item{fn}{File path to store stan model file.}

\item{parallel_chains}{Number of parallel chains to run. Defaults to 4.}

\item{iter_warmup}{Number of warmup iterations before sampling.
Defaults to 1,000.}

\item{iter_sampling}{Number of samples to draw from each chain.
Defaults to 1,000.}

\item{variables}{A character vector declaring which sampled variables
to return.
Defaults to "poi" which is the count of the
protein of interest.}
}
\value{
A stanfit object and a matrix of draws stored in a named list.
}
\description{
Writes a stan model and draws samples from the prior distribution of
protein count given a vector of isotype control (igg) count.
}
\examples{
\dontrun{
n <- 100L
igg <- rnbinom(100, mu=50, size=2)

prior_samples <- ADTnorm_prior_predict(n=n, igg=igg, variables="poi")
mu_poi <- apply(prior_samples$draws, 2, mean)
plot(igg, mu_poi)
}
}
\references{
\url{https://mc-stan.org/docs/reference-manual/}
}
