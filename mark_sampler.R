
# sampler
library(MCMCpack)
valc	<- read.table('valc.txt', header = TRUE)
t		<- valc$t
t <- readRDS("tweetstorm_data.rds") #; t <- (t*1)[1:400] # truncating data
t <- readRDS("tweet_storms.rds") %>% filter(!post_election) %>% pull(days_elapsed)


fcth	<- function(th, lam, t){
	n		<- length(t)
	ft		<- (th^n)/((lam^th)^(n+1)) # corrected n-1 to n+1, mark confirmed
	sti		<- (th - 1)*sum(log(t))
	st		<- exp(sti)
	tti		<- -sum((t/lam)^th)
	tt		<- exp(tti)
	out		<- ft*st*tt
	# out		<- (th^n)/((lam^th)^(n+1)) * exp((th - 1)*sum(log(t))) * exp(-sum((t/lam)^th)) # just for visibility
	return(out)
}

logfcth <- function(theta, lambda, data){
	n <- length(data)
	return(
		n*log(theta) - (n+1)*theta*log(lambda) + (theta-1)*sum(log(data)) - (1/(lambda^theta))*sum(data^theta)
	)
}

lam 		<- 100
th_grid		<- seq(0.001, 2, by = 0.001)
fcthplot	<- vector('numeric', length = length(th_grid))
for(i in 1:length(th_grid)){
	th			<- th_grid[i]
	fcthplot[i]	<- fcth(th, lam, t)
}

plot(th_grid, fcthplot, type = 'l')

# B		<- 50000
# a1		<- 4.8
# b1		<- 3.2
# # thet	<- 0.1

weibullSamp	<- function(t, B, a1, b1, thet, lam0){
	n		<- length(t)
	lambda	<- vector('numeric', length = B)
	theta	<- vector('numeric', length = B)
	ar		<- vector('numeric', length = B)
	lambda[1]	<- lam0
	theta[1]	<- thet
	
	for(b in 2:B){
		
		## sample lambda ##
		beta		<- sum(t^theta[b-1])
		lamth		<- rinvgamma(1, n, beta)
		lambda[b]	<- exp(log(lamth)*(theta[b-1]^(-1)))
		
		## sample theta ##
		tstar		<- rgamma(1, a1, b1)
		# aprob		<- min(1, (fcth(tstar, lambda[b-1], t)/fcth(thet, lambda[b-1], t))/(dgamma(tstar, a1, b1)/dgamma(thet, a1, b1)))
		aprob		<- min(1, 
									exp( logfcth(tstar, lambda[b-1], t) - logfcth(thet, lambda[b-1], t) ) / 
										(dgamma(tstar, a1, b1)/dgamma(thet, a1, b1))) # log space transform
		u			<- runif(1)
		if(u < aprob){
			thet	<- tstar
			ar[b]	<- 1
		}
		theta[b]	<- thet
	}
	out			<- NULL
	out$ar		<- ar[-(1:(B/2))]
	out$theta	<- theta[-(1:(B/2))]
	out$lambda	<- lambda[-(1:(B/2))]
	
	return(out)
}

lapply(weibullSamp(t, 10000, 3, 3, 0.5, 2), mean)

test <- weibullSamp(t, 100000, 136, 160, 0.5, 2)
mean(test$ar)
mcmcplot2(test$theta, "theta")
mcmcplot2(test$lambda, "lambda")


# visualizer (matrix shows x on row, y on column, in increasing order; lighter = greater)
# th_grid		<- seq(0.01, 3, by = 0.01) # mark original
# lam 		<- seq(7.01, 10, by = .01) # mark original
th_grid	<- seq(0.875, 0.975, length.out=500)
lam <- seq(1.650, 2.025, length.out=500)
fcthplot <- matrix(0, nrow = length(th_grid), ncol = length(lam))
for(i in 1:length(th_grid)){
	th <- th_grid[i]
	for(j in 1:length(lam)){
		lm	<- lam[j]
		fcthplot[i,j]	<- fcth(th, lm, t)
	}
}
image(fcthplot)

quantile(seq(0.75, 1.25, length.out=500), c(0.25, 0.45)) # 0.93
quantile(seq(1.5, 2.25, length.out=500), c(0.2, 0.7)) # 1.8

dev.off()

# Function to tune proposal density based on acceptance rate
# weibullSamp	<- function(t, B, a1, b1, thet, lam0)
tune_acceptance_rate <- function(a_vals, b_vals){
	acceptance_grid <- matrix(NA, nrow = length(a_vals), ncol = length(b_vals))
	for (i in 1:nrow(acceptance_grid)){
		for (j in 1:ncol(acceptance_grid)){
			acceptance_grid[i, j] <- mean(weibullSamp(t = t, 
																					 B = 10000,
																					 a1 = a_vals[i],
																					 b1 = b_vals[j],
																					 thet = 0.5,
																					 lam0 = 2)$ar)
		}
	}
	dimnames(acceptance_grid) <- list(a_vals, b_vals)
	return(acceptance_grid)
}

(tune_results <- tune_acceptance_rate(
	a_vals = seq(130, 140, by=2), 
	b_vals = seq(130, 140, by=2)/0.85)) # Choose Gamma(,)


a_vals = seq(2, 5, by=1) 
b_vals = seq(2, 5, by=1)
i <- 3; j <- 3

mean(weibullSamp(t = t, 
						B = 20000,
						a1 = 50,
						b1 = 50/0.85,
						thet = 0.5,
						lam0 = 2)$ar)
