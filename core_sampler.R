# Script for core Gibbs/MH sampler

# Load data
data <- readRDS("tweetstorm_data.rds")

# Sampler implementation
##  Full conditional theta
fcth <- function(lambda, theta, data) {
	n <- length(data)
	return(
		theta^n * lambda^(-theta*(n+1)) * 
			prod(data)^(theta-1) * 
			exp(-(1/lambda)^theta * sum(data^theta)))
}

## Implement sampler
execute_GibbsMH <- function(B, seed, start_lambda, start_theta, tune_a, tune_b){
	
	# DEBUG ONLY
	B = 10000 
	seed = 1234 
	start_lambda = 1.07 
	start_theta = 2 
	tune_a = 1.07 
	tune_b = 1
	
	set.seed(seed)
	
	lambda	<- vector("numeric", B)
	theta <- vector("numeric", B)
	ar <- vector("numeric", B)
	
	lambda[1] <- start_lambda
	th <- start_theta
	theta[1] <- th
	
	n	<- length(data)
	
	for(t in 2:B){
		# DEBUG ONLY
		t <- 2
		
		### sample lambda ###
		lambda[t] <- rgamma(1, n+2, sum(data^theta[t-1]))^(-1/theta[t-1])
		
		### sample theta ###
		thstar <- rgamma(1, tune_a, tune_b)
		aprob  <- min(1, 
									( fcth(lambda[t-1], thstar, data) / 
											fcth(lambda[t-1], theta[t-1], data) ) /
										( dgamma(thstar, tune_a, tune_b) / 
												dgamma(theta[t-1], tune_a, tune_b) ) )
		U <- runif(1)
		if(U < aprob){
			th <- thstar
			ar[t] <- 1
		}
		theta[t] <- th
	}
	
	return(list(
		lambda = lambda[-c(1:(B/2))], 
		theta = theta[-c(1:(B/2))], 
		ar = mean(ar)))
}

execute_GibbsMH(
	B = 10000,
	seed = 1234,
	start_lambda = 1.07,
	start_theta = 2,
	tune_a = 1.07,
	tune_b = 1)
# RESULT: Error in if (U < aprob) { : missing value where TRUE/FALSE needed
# Reason: fcth (full conditional theta) = NaN