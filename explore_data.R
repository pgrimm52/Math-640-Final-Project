# Preliminary exploration of dataset

# load libraries
library(jsonlite)
library(dplyr)
library(lubridate)
library(stringr)
library(ggplot2)
library(dbscan)

# load data
data <- readRDS("trump_tweets.rds") %>%
	tbl_df() %>%
	arrange(desc(date)) %>%
	filter(!is_retweet, 
				 source %in% c("Twitter for iPhone", "Twitter for Android"),
				 year(date) > 2015) %>%
	mutate(previous = lead(date),
				 hours_elapsed = (date - previous)/3600) %>%
	slice(-6285) # omit last entry with blank hours_elapsed
	
# seconds -> hours: div by 3600
# seconds -> days: div by 86400

summary(data$hours_elapsed)
quantile(data$hours_elapsed, c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975), na.rm=TRUE)

qplot(hours_elapsed, data=data, facets = year(date) ~ .)

summary(as.numeric(data$hours_elapsed)[year(data$date) > 2017])

# explore data
data %>% 
	count(source)

qplot(date, 
			data = data %>% 
				filter(source %in% c("Twitter Web Client", "Twitter for iPhone", "Twitter for Android")), 
			geom="histogram", 
			facets = source ~ . )

data %>%
	count(is_retweet)

# clustering tweetstorms
test_subset <- data %>% slice(1:1000)
test_subset <- data

test_subset$cluster <- dbscan(as.matrix(test_subset$date), eps = 1000, minPts = 3)$cluster

test_subset %>%
	count(cluster)

pilot_data <- test_subset %>%
	filter(cluster != 0) %>%
	group_by(cluster) %>%
	summarize(date = mean(date)) %>%
	mutate(previous = lead(date),
				 time_elapsed = as.duration(date - previous)/ddays(1)) # time_elapsed converted to days

qplot(time_elapsed, data=pilot_data)

##################
# Weibull sampler
##################

# Load data
data <- pilot_data$time_elapsed[!is.na(pilot_data$time_elapsed)]

# Full conditional theta
fcth <- function(lambda, theta, data) {
	n <- length(data)
	return(
		theta^n * lambda^(-theta*(n+1)) * 
			prod(data)^(theta-1) * 
			exp(-(1/lambda)^theta * sum(data^theta)))
}

# Implement sampler
execute_GibbsMH_q2 <- function(B, seed, start_lambda, start_theta, tune_a, tune_b){
	
	set.seed(seed)
	
	lambda	<- vector("numeric", B)
	theta <- vector("numeric", B)
	ar <- vector("numeric", B)
	
	lambda[1] <- start_lambda
	th <- start_theta
	theta[1] <- th
	
	n	<- length(data)
	
	for(t in 2:B){
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

tune_acceptance_rate <- function(a_vals, b_vals){
	
	acceptance_grid <- matrix(NA, nrow = length(a_vals), ncol = length(b_vals))
	for (i in 1:nrow(acceptance_grid)){
		for (j in 1:ncol(acceptance_grid)){
			acceptance_grid[i, j] <- execute_GibbsMH_q2(B = 50000, seed = 121, 
																									start_lambda = 1, 
																									start_theta = 0.1, 
																									tune_a = a_vals[i], 
																									tune_b = b_vals[j])$ar
		}
	}
	dimnames(acceptance_grid) <- list(a_vals, b_vals)
	return(acceptance_grid)
}

# Proposal distribution
hist(sampler_data)
plot(sampler_data)
curve(dgamma(x, 5, 2), col="blue", lty=2, lwd=2, add=TRUE)
curve(dgamma(x, 200, 1), col="blue", lty=2, lwd=2, xlim=c(0,500))


# Execute sampler
chain1 <- execute_GibbsMH_q2(
	B = 10000, seed = 1234, 
	start_lambda = 1, start_theta = 0.1, 
	tune_a = 200, tune_b = 1)
chain1$ar
