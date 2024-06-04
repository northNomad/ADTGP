set.seed(123)
n <- 100
rho <- .3
runif(n, 0, 1)


mu_poi <- 500
thi_poi <- 5

mu_igg <- 20
thi_igg <- 2

library(cmdstanr)
set_cmdstan_path()

MASS::mvrnorm(n, mu = c(0, 0),
              Sigma = matrix(c(1, rho, rho, 1), nrow=2, byrow = FALSE)) -> d

cor(d[, 1], d[, 2]) #.3295
plot(d[, 1], d[, 2])

d_cdf <- list(d1=ecdf(d[, 1]), d2=ecdf(d[, 2]))

d1_cd <- sapply(d[, 1], function(x) d_cdf$d1(x))
d2_cd <- sapply(d[, 2], function(x) d_cdf$d2(x))
plot(d1_cd, d2_cd)
cor(d1_cd, d2_cd) #.3227

dist_poi <- rnbinom(1e5, mu=mu_poi, size=thi_poi)
dist_igg <- rnbinom(1e5, mu=mu_igg, size=thi_igg)
plot(dist_igg, dist_poi)

n_poi <- quantile(dist_poi, probs = d1_cd)
n_igg <- quantile(dist_igg, probs = d2_cd)

plot(n_igg, n_poi)


m <- cmdstanr::cmdstan_model("fit_model.stan")
data_sim <- list(N=n, igg=n_igg, poi=n_poi, nrep=100, igg_rep=rep(0, 100))
