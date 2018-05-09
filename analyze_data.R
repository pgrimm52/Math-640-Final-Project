########################
# File: analyze_data.R
# Purpose: to conduct core sampling and analysis
# Output: console, graphs
########################

###############
# A. Setup 
###############

rm(list=ls())

# Load packages and helper functions
source("helper_functions.R")

# Load data
tweet_storms  <- readRDS("tweet_storms.rds")

all_storms    <- tweet_storms %>% pull(days_elapsed)
pre_election  <- tweet_storms %>% filter(!post_election) %>% pull(days_elapsed)
post_election <- tweet_storms %>% filter(post_election) %>% pull(days_elapsed)

###############
# B. Examine different posteriors 
# (always source sampler file before runs, since posteriors overwrite each other) 
###############

# LOG-NORMAL #############################
source("lognormal_sampler.R")

# NON-INFORMATIVE
# 1. No tuning necessary
# 2. Sample 
lognorm_samples_noninf <- lognormSamp(data=all_storms, B=10000, 
															 a=1000, b=1000, 
															 c=1, d=1000)

# 3. Check for convergence
mcmcplot2(lognorm_samples_noninf$mu)
mcmcplot2(lognorm_samples_noninf$sig2)

# 4. Grab Log-normal parameters
lapply(lognorm_samples_noninf, mean)
lapply(lognorm_samples_noninf, median)

# INFORMATIVE
# 1. No tuning necessary
# 2. Sample 
lognorm_samples_inf <- lognormSamp(data=all_storms, B=10000, 
															 a=2, b=1, 
															 c=0.752763, d=1)

# 3. Check for convergence
mcmcplot2(lognorm_samples_inf$mu)
mcmcplot2(lognorm_samples_inf$sig2)

# 4. Grab Log-normal parameters
lapply(lognorm_samples_inf, mean)
lapply(lognorm_samples_inf, median)


# GAMMA #############################
source("gamma_sampler.R")

# NON-INFORMATIVE
# 1. Tune with Gamma proposal
tune_acceptance_rate(
	a_vals = seq(60, 70, by=1),
	b_vals = seq(60, 70, by=1),
	sampler = gammaSamp,
	data = all_storms,
	B = 10000,
	p=1, q=0, r=0, s=0,
	alpha_start = 1, beta_start = 1) # settle on Gamma(60, 63) with ~43% acceptance

# 2. Sample 
gamma_samples_noninf <- gammaSamp(data = all_storms, B = 10000,
													 a1 = 60, b1 = 63,
													 p=1, q=0, r=0, s=0,
													 alpha_start = 1, beta_start = 1)

# 3. Check for convergence
mcmcplot2(gamma_samples_noninf$alpha)
mcmcplot2(gamma_samples_noninf$beta)

# 4. Grab Gamma parameters
lapply(gamma_samples_noninf, mean)
lapply(gamma_samples_noninf, median)

# INFORMATIVE
# 1. Tune with Gamma proposal
# Not necessary

# 2. Sample 
gamma_samples_inf <- gammaSamp(data = all_storms, B = 10000,
													 a1 = 64, b1 = 67,
													 p=10, q=10, r=10, s=10,
													 alpha_start = 1, beta_start = 1)

# 3. Check for convergence
mcmcplot2(gamma_samples_inf$alpha)
mcmcplot2(gamma_sample_infs$beta)

# 4. Grab Gamma parameters
lapply(gamma_samples_inf, mean)
lapply(gamma_samples_inf, median)


# WEIBULL #############################
source("weibull_sampler.R")

# 1. Tune with Gamma proposal
tune_acceptance_rate(
	a_vals = seq(60, 70, by=1),
	b_vals = seq(60, 70, by=1),
	sampler = weibullSamp,
	data = all_storms,
	B = 20000,
	theta_start = 1, lambda_start = 1) # settle on Gamma(60, 67) with ~44% acceptance

# 2. Sample 
weibull_samples <- weibullSamp(
	seed = 1, data = all_storms, B = 20000,
	a1 = 60, b1 = 67,
	theta_start = 0.5, lambda_start = 1)

# 3. Check for convergence
mcmcplot2(weibull_samples$theta)
mcmcplot2(weibull_samples$lambda)

# 4. Grab Weibull parameters
lapply(weibull_samples, mean)
lapply(weibull_samples, median)

# OPTIONAL EXHAUSTIVE CONVERGENCE CHECK CODE
# # X. Run multiple chains and check for convergence
# chain1 <- weibullSamp(seed = 1, data = all_storms, B = 20000,
# 											a1 = 60, b1 = 67,
# 											theta_start = 0.5, lambda_start = 1)
# chain2 <- weibullSamp(seed = 2, data = all_storms, B = 20000,
# 											a1 = 60, b1 = 67,
# 											theta_start = 0.75, lambda_start = 2)
# chain3 <- weibullSamp(seed = 3, data = all_storms, B = 20000,
# 											a1 = 60, b1 = 67,
# 											theta_start = 1, lambda_start = 3)
# chain4 <- weibullSamp(seed = 4, data = all_storms, B = 20000,
# 											a1 = 60, b1 = 67,
# 											theta_start = 1.5, lambda_start = 4)
# 
# gelman_rubin(chain1$theta, chain2$theta, chain3$theta, chain4$theta)
# gelman_rubin(chain1$lambda, chain2$lambda, chain3$lambda, chain4$lambda)
# 
# # If desired: convergence plots for every chain variable
# # mcmcplot2(chain1$theta) # etc.
# 
# # Y. Determine optimal thinning (10 seems best)
# show_thinning_options(chain1$theta)
# show_thinning_options(chain2$theta)
# show_thinning_options(chain3$theta)
# show_thinning_options(chain4$theta)
# 
# show_thinning_options(chain1$lambda)
# show_thinning_options(chain2$lambda)
# show_thinning_options(chain3$lambda)
# show_thinning_options(chain4$lambda)


###############
# C. Decide on best likelihood (probably tie between weibull, gamma)
###############

# Graphical comparison to empirical density
plot(density(all_storms, from=0), lwd=2, ylim=c(0, 0.4))
curve(dlnorm(x, 0.4351, sqrt(1.086)), col="green", lwd=2, lty=2, add=TRUE)
curve(dgamma(x, 0.9197, 0.3163), col="blue", lwd=2, lty=2, add=TRUE)
curve(dweibull(x, 0.9083, 2.8735), col="red", lwd=2, lty=2, add=TRUE)

# Show replicated moment distribution with original data
show_replicate_analysis(all_storms, rlnorm, lognorm_samples$mu, sqrt(lognorm_samples$sig2)) # remember sqrt!
show_replicate_analysis(all_storms, rgamma, gamma_samples$alpha, gamma_samples$beta)
show_replicate_analysis(all_storms, rweibull, weibull_samples$theta, weibull_samples$lambda)
dev.off()

###############
# D. Compare credible intervals of pre- and post-election 
# (there is a difference, but not stat significant because intervals overlap)
###############

# Use Weibull because paramters have interesting interpretation
source("weibull_sampler.R")

weibull_samples_pre <- weibullSamp(
	seed = 1, data = pre_election, B = 20000,
	a1 = 60, b1 = 67,
	theta_start = 0.5, lambda_start = 1)

weibull_samples_post <- weibullSamp(
	seed = 1, data = post_election, B = 20000,
	a1 = 60, b1 = 67,
	theta_start = 0.5, lambda_start = 1)

posterior_summary(weibull_samples_pre$theta)
posterior_summary(weibull_samples_post$theta)

posterior_summary(weibull_samples_pre$lambda)
posterior_summary(weibull_samples_post$lambda)

# Look at median of distribution via parameters
posterior_summary(weibull_samples_pre$lambda * log(2)^(1/weibull_samples_pre$theta))
posterior_summary(weibull_samples_post$lambda * log(2)^(1/weibull_samples_post$theta))

###############
# E. How well does pre-election model/data do for post-election realized data?
###############

weibull_samples_pre <- weibullSamp(
	seed = 1, data = pre_election, B = 20000,
	a1 = 60, b1 = 67,
	theta_start = 0.5, lambda_start = 1)

posterior_summary(weibull_samples_pre$theta)
posterior_summary(weibull_samples_pre$lambda)

plot(density(post_election, from=0), lwd=2, ylim=c(0, 0.4))
curve(dweibull(x, 0.8378, 3.4118), add=TRUE, col="red", lty=2, lwd=2)
