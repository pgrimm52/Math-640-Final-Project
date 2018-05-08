#######################
# Gamma sampler code
#######################

#### CONJUGATE prior

# Posterior
post <- function(data, p, q, r, s, alpha, beta){
	n <- length(data)
	return(
		beta^(alpha*(s+n)) / gamma(alpha)^(r+n) * (p * prod(data))^(alpha-1)
	)
}

# Log-Posterior
logpost <- function(data, p, q, r, s, alpha, beta){
	n <- length(data)
	return(
		(alpha*(s+n))*log(beta) - (r+n)*log(gamma(alpha)) + (alpha-1)*(log(p) + sum(log(data)))
	)
}

# Sampler implementation
# to test: gammaSamp(data=data, B=10000, p=1, q=1, r=1, s=1, a1=1, b1=1, alpha_start=1, beta_start=1)
gammaSamp	<- function(seed=123, data, B, 
											p, q, r, s, 
											a1, b1, 
											alpha_start, beta_start){
	
	set.seed(seed)
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
									exp( logpost(data=data, p=p, q=q, r=r, s=s, alpha = astar, beta = beta[b-1]) -
										logpost(data=data, p=p, q=q, r=r, s=s, alpha = alpha[b-1], beta = beta[b-1]) ) /
									( dgamma(astar, a1, b1) /
										dgamma(alpha[b-1], a1, b1) ) )
		u	<- runif(1)
		alpha[b]	<- alpha[b-1]
		if(u < aprob){
			alpha[b]	<- astar
			ar[b]	<- 1
		}
	}

	return(list(
		ar = ar[-(1:(B/2))],
		alpha	= alpha[-(1:(B/2))],
		beta = beta[-(1:(B/2))]
	))
}
