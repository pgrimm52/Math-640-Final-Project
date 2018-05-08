######################################################
# HELPER FUNCTIONS with applicablility across scripts
# source("helper_functions.R")
######################################################

#####################################
# 0. IMPORTANT PACKAGES
library(jsonlite)
library(dplyr)
library(lubridate)
library(stringr)
library(ggplot2)
library(dbscan)
library(mcmcplots)
library(MCMCpack)
library(moments)

# hist(all_storms)
# curve(dwhatever(x, a, b), col="blue", add=TRUE)

#####################################
# 1. ACCEPTANCE RATE TUNER
# (for Gamma proposals only; requires sampler with specific style)
tune_acceptance_rate <- function(a_vals, b_vals, sampler, ...){
	acceptance_grid <- matrix(NA, nrow = length(a_vals), ncol = length(b_vals))
	for (i in 1:nrow(acceptance_grid)){
		for (j in 1:ncol(acceptance_grid)){
			acceptance_grid[i, j] <- mean(sampler(a1 = a_vals[i], b1 = b_vals[j], ...)$ar)
		}
	}
	dimnames(acceptance_grid) <- list(a_vals, b_vals)
	return(acceptance_grid)
}

# EXAMPLE USAGE
# tune_acceptance_rate(
# 	a_vals = seq(40, 42, by=1), 
# 	b_vals = seq(40, 42, by=1),
# 	sampler = gammaSamp,
# 	data = t, 
# 	B = 10000,
# 	p = 1, q = 1, r = 1, s = 1, 
# 	alpha_start = 1, 
# 	beta_start = 1
# )

#####################################
# 2. MCMCPLOT WRAPPER
mcmcplot2 <- function(samp, name="variable"){
	mat <- matrix(samp, ncol=1)
	colnames(mat) <- name
	mcmcplots::mcmcplot1(mat)
}

#####################################
# 3. WEIBULL SUMMARIZER
weibull_summary <- function(k, lam){
	mean <- lam * gamma(1 + 1/k)
	median <- lam * log(2)^(1/k)
	variance <- lam^2 * (gamma(1 + 2/k) - (gamma(1 + 1/k))^2)
	return(list(mean = mean, 
							median = median, 
							variance = variance))
}

#####################################
# 4. POSTERIOR SUMMARIZER
posterior_summary <- function(data){
	return(list(
		mean = mean(data),
		quantiles = quantile(data, c(0.025, 0.5, 0.975))
	))
}

#####################################
# 5. GELMAN-RUBIN TESTER
gelman_rubin <- function(chain1, chain2, chain3, chain4){
	allChains <- mcmc.list(list(
		mcmc(chain1), 
		mcmc(chain2), 
		mcmc(chain3),
		mcmc(chain4)))
	return(gelman.diag(allChains))
}

#####################################
# 6. THINNING TOOLKIT
thin <- function(data, x) {
	data[seq(x, length(data), by=x)]
}

show_thinning_options <- function(data){
	par(mfrow=c(1,4))
		acf(data)
		acf(thin(data, 5))
		acf(thin(data, 10))
		acf(thin(data, 20))
}

#####################################
# 7. REPLICATED DATA ANALYZER
show_replicate_analysis <- function(data, rdistribution, sampled_parameter1, sampled_parameter2){
	
	replicated_data <- apply(
		cbind(sampled_parameter1, sampled_parameter2),
		1,
		function(x){ rdistribution(length(data), x[1], x[2]) })
	
	replicated_stats <- apply(
		replicated_data,
		2,
		function(x) { c(median(x), mean(x), var(x), skewness(x), kurtosis(x)) })
	
	original_stats <- do.call(
		function(x) { c(median(x), mean(x), var(x), skewness(x), kurtosis(x)) },
		list(data))
	
	pval <- sapply(1:5, function(i){
		sum(replicated_stats[i, ] < original_stats[i])/length(replicated_stats[i, ])
	})
	
	par(mfrow=c(2,3))
	hist(replicated_stats[1, ], main = "Median"); abline(v = original_stats[1], col = "red"); legend("topright", paste("P-value:\n", pval[1], sep=""), bty="n") 
	hist(replicated_stats[2, ], main = "Mean"); abline(v = original_stats[2], col = "red"); legend("topright", paste("P-value:\n", pval[2], sep=""), bty="n")
	hist(replicated_stats[3, ], main = "Variance"); abline(v = original_stats[3], col = "red"); legend("topright", paste("P-value:\n", pval[3], sep=""), bty="n")
	hist(replicated_stats[4, ], main = "Skewness"); abline(v = original_stats[4], col = "red"); legend("topright", paste("P-value:\n", pval[4], sep=""), bty="n")
	hist(replicated_stats[5, ], main = "Kurtosis"); abline(v = original_stats[5], col = "red"); legend("topright", paste("P-value:\n", pval[5], sep=""), bty="n")
}