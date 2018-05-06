tweetstorm <- readRDS("/Users/ibrlouise/Documents/Georgetown Grad School/Grad School - Spring 2018/Bayesian Statistics/tweet_storms.rds")
tt <- tweetstorm$days_elapsed

library(MCMCpack)
library(mcmcplots)

log_normal_gibbs <- function(B, seed, data, a, b, c, d){
  n <- length(data)
  mu <- mean(data)
  sig <- var(data)
  
  sigmas <- matrix(NA, nrow=B, ncol=1)
  mus <- matrix(NA, nrow=B, ncol=1)
  
  mus[1,] <- mu
  sigmas[1,] <- sig
  
  set.seed(seed)
  
  for(t in 2:B){
    mus[t, 1] <- rnorm(1, 
                       (((sum(log(data)))/sig) + (c/d))
                       /((n/sig)+(1/d)), 
                       ((n/sig)+(1/d))^(-1))
    sigmas[t, 1] <- rinvgamma(1, (n/2)+a, (1/2)*(sum((log(data)-mu)^2)))
  }
  
  mumean <- mean(mus[-(1:(B/2))])
  sigmean <- mean(sigmas[-(1:(B/2))])
  
  mupPlot	<- mcmc.list(list(mcmc(mus)))
  sigPlot <- mcmc.list(list(mcmc(sigmas)))
  
  varnames(mupPlot) <- c("mu")
  varnames(sigPlot) <- c("sigma^2")
  
  return(c("mean mu"= mumean, "mean sigma" = sigmean, mcmcplot1(mupPlot, greek=TRUE), mcmcplot1(sigPlot, greek=TRUE)))
  
}

log_normal_gibbs(
  B = 10000,
  seed = 1234,
  data = tt,
  a = 1,
  b = 1,
  c = 2,
  d = 1
)


plot(density(tt, from=0))
curve(dlnorm(x, 0.5101603, 8.3773659), add=TRUE, col="blue")