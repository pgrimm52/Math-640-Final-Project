# Gamma sampler script

# Load data
t <- readRDS("tweet_storms.rds") %>% pull(days_elapsed)

# Gibbs/M-H implementation
post <- function(data, p, q, r, s, alpha, beta){
	n <- length(data)
	return(
		beta^(alpha*(s+n)) / gamma(alpha)^(r+n) * (p * prod(data))^(alpha-1)
	)
}

logpost <- function(data, p, q, r, s, alpha, beta){
	n <- length(data)
	return(
		(alpha*(s+n))*log(beta) - (r+n)*log(gamma(alpha)) + (alpha-1)*(log(p) + sum(log(data)))
	)
}

gammaSamp	<- function(data, B, p, q, r, s, a1, b1, alpha_start, beta_start){
	n	<- length(data)
	alpha	<- vector('numeric', length = B)
	beta	<- vector('numeric', length = B)
	ar <- vector('numeric', length = B)
	alpha[1]	<- alpha_start
	beta[1]	<- beta_start
	
	for(b in 2:B){
		
		## sample beta ##
		beta[b] <- rgamma(1, alpha[b-1]*(s + n) + 1, q + sum(data))	

		## sample alpha ##
		astar		<- rgamma(1, a1, b1)
		# aprob		<- min(1, 
		# 							( post(data, p, q, r, s, alpha = astar, beta = beta[b-1]) / 
		# 								post(data, p, q, r, s, alpha = alpha[b-1], beta = beta[b-1]) ) /
		# 							( dgamma(astar, a1, b1) /
		# 								dgamma(alpha[b-1], a1, b1) ) )
		aprob		<- min(1,
									exp( logpost(data, p, q, r, s, alpha = astar, beta = beta[b-1]) -
										logpost(data, p, q, r, s, alpha = alpha[b-1], beta = beta[b-1]) ) /
									( dgamma(astar, a1, b1) /
										dgamma(alpha[b-1], a1, b1) ) )
		u	<- runif(1)
		alpha[b]	<- alpha[b-1]
		if(u < aprob){
			alpha[b]	<- astar
			ar[b]	<- 1
		}
	}

	out				<- NULL
	out$ar		<- ar[-(1:(B/2))]
	out$alpha	<- alpha[-(1:(B/2))]
	out$beta	<- beta[-(1:(B/2))]
	return(out)
}

result <- gammaSamp(
		data = data, 
		B = 10000, 
		p = 1, q = 1, r = 1, s = 1, 
		a1 = 5, b1 = 1, 
		alpha_start = 1, 
		beta_start = 1)

lapply(result, mean)

hist(t)
curve(5e2*dgamma(x, 1.04, 0.415), add=TRUE, col="red")

# Tuning

tune_acceptance_rate <- function(a_vals, b_vals){
	acceptance_grid <- matrix(NA, nrow = length(a_vals), ncol = length(b_vals))
	for (i in 1:nrow(acceptance_grid)){
		for (j in 1:ncol(acceptance_grid)){
			acceptance_grid[i, j] <- mean(gammaSamp(data = t, 
																							B = 10000,
																							p = 1, q = 1, r = 1, s = 1, 
																							a1 = a_vals[i], b1 = b_vals[j], 
																							alpha_start = 1, 
																							beta_start = 1)$ar)
		}
	}
	dimnames(acceptance_grid) <- list(a_vals, b_vals)
	return(acceptance_grid)
}

(tune_results <- tune_acceptance_rate(
	a_vals = seq(30, 40, by=1), 
	b_vals = seq(30, 40, by=1)))

result <- gammaSamp(
	data = data, 
	B = 40000, 
	p = 1, q = 1, r = 1, s = 1, 
	a1 = 34, b1 = 38, 
	alpha_start = 1, 
	beta_start = 1)
