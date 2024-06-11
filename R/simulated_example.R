#' Runs the simulated example in the manuscript
#'
#' Returns the true dataframe, the noisy dataframe, and the stanfit object.
#' Also draws two plots in the working directory as a side effect.
ADTGP_RunSimulation <- function(){
  set.seed(6)
  ## parameters for poi count
  mu1_poi <- 200
  mu2_poi <- 100
  phi_poi <- 1

  ## parameters for igg count
  mu_igg <- 20
  phi_igg <- 5

  ## number of cells to simulate
  n <- 50

  ## create distribution
  big_poi1 <- rnbinom(1e4, mu=mu1_poi, size=phi_poi)
  big_poi2 <- rnbinom(1e4, mu=mu2_poi, size=phi_poi)
  big_igg <- rnbinom(1e4, mu=mu_igg, size=phi_igg)

  ## Draw correlated count data
  sim_quantile1 <- runif(n, 0, 1)
  sim_quantile2 <- runif(n, 0, 1)

  n1_poi <-quantile(big_poi1, prob=sim_quantile1)
  n1_igg <- quantile(big_igg, prob=sim_quantile1)
  n2_poi <- quantile(big_poi2, prob=sim_quantile2)
  n2_igg <- quantile(big_igg, prob=sim_quantile2)

  ## Put data in dataframe 'd'
  data.frame(
    T = rep(c(1, 2), each=n),
    poi = c(n1_poi, n2_poi ),
    igg = c(n1_igg, n2_igg)
  ) -> d

  ## Sample from data.
  ## For some reason, we:
  ## Observe more low igg count cells in T=1
  ## Observe more high igg count cells in T=2
  ## This masks true variation between T=1 and T=2
  index1 <- sample(1:length(n1_poi), size=50, replace = TRUE, prob=1/n1_igg)
  index2 <- sample(1:length(n2_poi), size=50, replace = TRUE, prob=n2_igg)

  ###########
  ## Plot: Density plot of protein X expression:
  ###########
  svg("p_ProteinXExpressionDensity.svg", width=5.6, height=4.30)
  plot(NULL, xlim=c(0, 9), ylim=c(0, .6), ylab="Density", xlab="log (n+1)",
       main="Protein X Expression")
  lines(density(log(n1_poi + 1)), col="red", lty=1)
  lines(density(log(n1_poi[index1] + 1)), col="red", lty=2)
  lines(density(log(n2_poi + 1)), col="blue", lty=1)
  lines(density(log(n2_poi[index2] + 1)), col="blue", lty=2)
  legend(x=0, y=.55,
         legend=c("Truth (T=1)", "Observed (T=1)", "Truth (T=2)", "Observed (T=2)"),
         col=c("red", "red", "blue", "blue"),
         lty=c(1, 2, 1, 2),
         cex=.8)
  dev.off()
  ###########

  ## Put observed data in data.frame 'dsub'
  dsub <- data.frame(
    T=rep(c(1, 2), times=c(length(index1), length(index2))),
    poi=c(n1_poi[index1], n2_poi[index2]),
    igg=c(n1_igg[index1], n2_igg[index2])
  )

  ## Prepare design matrix and fit model with ADTGP
  dm <- matrix(nrow=nrow(dsub), ncol=2)
  colnames(dm) <- c("mu0", "T")
  dm[, "mu0"] <- 1
  dm[, "T"] <- dsub$T

  fit <- ADTGP(igg=dsub$igg, poi=dsub$poi, design_matrix=dm,
              iter_warmup = 1500, iter_sampling = 1500)

  ## Get posterior of Treatment coefficient and calculate l2fc
  bT_post <- fit$stanfit$draws("b_T", format="draws_matrix")
  post_l2fc <- log2(exp(bT_post[, "b_T[1]"] - bT_post[, "b_T[2]"]))

  ## Calculate observed l2fc
  obs_mu1 <- mean(dsub[dsub$T == 1, "poi"])
  obs_mu2 <- mean(dsub[dsub$T == 2, "poi"])
  obs_l2fc <- log2(obs_mu1 / obs_mu2)

  ## Calculate true l2fc
  truth_l2fc <- log2(mean(n1_poi) / mean(n2_poi))

  ###########
  ## Plot: Boxplot of L2FC between T=1 and T=2
  ###########
  svg("p_boxplotL2FC.svg", width=4, height=4.30)
  boxplot(post_l2fc, ylim=c(-.5, 1.4), ylab="Log2-Fold-Change",
          main="T=1 v. T=2 L2FC Posterior")
  points(x=1, y=obs_l2fc, col="red", pch=19, cex=1.3)
  points(x=1, y=truth_l2fc, col="blue", pch=19, cex=1.3)
  legend(x=.5, y=0,
         pch=c(19, 19),
         col=c("blue", "red"),
         legend=c("Truth", "Observed"))
  dev.off()
  ###########

  list(df_true = d,
       df_noisy = dsub,
       stanfit = fit) -> out
  return(out)
}
