#######################
# Log-Normal sampler code
#######################

#### NORMAL/INVGAMMA priors

# Sampler implementation
# to test: lognormSamp(data = data, B=10000, a=1, b=1, c=1, d=1)
lognormSamp <- function(seed=123, data, B, a, b, c, d){
  
	set.seed(seed)
	n <- length(data)

  sig2 <- vector(length = B)
  mu <- vector(length = B)
  
  mu[1] <- mean(data)
  sig2[1] <- var(data)
  
  for(i in 2:B){
    mu[i] <- rnorm(1,
    							 ( sum(log(data))/sig2[i-1] + c/d ) / ( n/sig2[i-1] + 1/d ),
    							 ( n/sig2[i-1] + 1/d )^(-1))
    sig2[i] <- rinvgamma(1, 
    										 a + n/2, 
    										 b + (1/2)*sum((log(data)-mu[i-1])^2))
  }
  
  return(list(
  	mu = mu[-(1:(B/2))],
  	sig2 = sig2[-(1:(B/2))]
  ))
}
