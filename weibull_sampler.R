#######################
# Weibull sampler code
#######################

#### NON-INFORMATIVE prior

# Posterior
post	<- function(th, lam, data){
	n		<- length(data)
	ft		<- (th^n)/((lam^th)^(n+1))
	sti		<- (th - 1)*sum(log(data))
	st		<- exp(sti)
	tti		<- -sum((data/lam)^th)
	tt		<- exp(tti)
	out		<- ft*st*tt # equivalent to (th^n)/((lam^th)^(n+1)) * exp((th - 1)*sum(log(data))) * exp(-sum((data/lam)^th))
	return(out)
}

# Log-Posterior
logpost <- function(theta, lambda, data){
	n <- length(data)
	return(
		n*log(theta) - (n+1)*theta*log(lambda) + (theta-1)*sum(log(data)) - (1/(lambda^theta))*sum(data^theta)
	)
}

# Sampler implementation
# to test: weibullSamp(data=data, B=10000, a1=1, b1=1, theta_start=1, lambda_start=1)
weibullSamp	<- function(seed=123, data, B, a1, b1, theta_start, lambda_start){
	
	set.seed(seed)
	n		<- length(data)
	lambda	<- vector('numeric', length = B)
	theta	<- vector('numeric', length = B)
	ar		<- vector('numeric', length = B)
	lambda[1]	<- lambda_start
	theta[1]	<- theta_start
	
	for(b in 2:B){
		
		## sample lambda ##
		beta		<- sum(data^theta[b-1])
		lamth		<- rinvgamma(1, n, beta)
		lambda[b]	<- exp(log(lamth)*(theta[b-1]^(-1)))
		
		## sample theta ##
		tstar		<- rgamma(1, a1, b1)
		# aprob		<- min(1, (post(tstar, lambda[b-1], data)/post(theta[b-1], lambda[b-1], data))/(dgamma(tstar, a1, b1)/dgamma(thet, a1, b1))) # 
		aprob		<- min(1, 
									exp(logpost(tstar, lambda[b-1], data) - logpost(theta[b-1], lambda[b-1], data)) / (dgamma(tstar, a1, b1)/dgamma(theta[b-1], a1, b1))) # log space transform
		u			<- runif(1)
		theta[b] <- theta[b-1]
		if(u < aprob){
			theta[b]	<- tstar
			ar[b]	<- 1
		}
	}

	return(list(
		ar = ar[-(1:(B/2))],
		theta	= theta[-(1:(B/2))],
		lambda = lambda[-(1:(B/2))]
	))
}

#### TO BE IMPLEMENTED: 
#### INFORMATIVE prior

# # Posteriors
# post_inf_theta <- function()
# post_inf_lambda <- function()