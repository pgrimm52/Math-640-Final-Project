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
															 a=1, b=10, 
															 c=1, d=1000)

# 3. Check for convergence
mcmcplot2(lognorm_samples_noninf$mu)
mcmcplot2(lognorm_samples_noninf$sig2)

# 4. Grab Log-normal parameters
posterior_summary(lognorm_samples_noninf$mu)
posterior_summary(lognorm_samples_noninf$sig2)

final_lognorm <- lognorm_samples_noninf

# INFORMATIVE
# 1. No tuning necessary
# 2. Sample 
lognorm_samples_inf <- lognormSamp(data=all_storms, B=10000, 
															 a=1, b=10, 
															 c=0.752763, d=1)

# 3. Check for convergence
mcmcplot2(lognorm_samples_inf$mu)
mcmcplot2(lognorm_samples_inf$sig2)

# 4. Grab Log-normal parameters
posterior_summary(lognorm_samples_inf$mu)
posterior_summary(lognorm_samples_inf$sig2)


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
gamma_samples_noninf.chain1 <- gammaSamp(
	seed = 1, data = all_storms, B = 20000,
	a1 = 60, b1 = 63,
	p=1, q=0, r=0, s=0,
	alpha_start = 0.7, beta_start = 0.2)

gamma_samples_noninf.chain2 <- gammaSamp(
	seed = 2, data = all_storms, B = 20000,
	a1 = 60, b1 = 63,
	p=1, q=0, r=0, s=0,
	alpha_start = 1.1, beta_start = 0.3)

gamma_samples_noninf.chain3 <- gammaSamp(
	seed = 3, data = all_storms, B = 20000,
	a1 = 60, b1 = 63,
	p=1, q=0, r=0, s=0,
	alpha_start = 0.9, beta_start = 0.4)

gamma_samples_noninf.chain4 <- gammaSamp(
	seed = 4, data = all_storms, B = 20000,
	a1 = 60, b1 = 63,
	p=1, q=0, r=0, s=0,
	alpha_start = 0.8, beta_start = 0.5)

# 3. Check for convergence
mcmcplot2(gamma_samples_noninf.chain1$alpha)
mcmcplot2(gamma_samples_noninf.chain2$alpha)
mcmcplot2(gamma_samples_noninf.chain3$alpha)
mcmcplot2(gamma_samples_noninf.chain4$alpha)

gelman_rubin(
	gamma_samples_noninf.chain1$alpha,
	gamma_samples_noninf.chain2$alpha,
	gamma_samples_noninf.chain3$alpha,
	gamma_samples_noninf.chain4$alpha)

mcmcplot2(gamma_samples_noninf.chain1$beta)
mcmcplot2(gamma_samples_noninf.chain2$beta)
mcmcplot2(gamma_samples_noninf.chain3$beta)
mcmcplot2(gamma_samples_noninf.chain4$beta)

gelman_rubin(
	gamma_samples_noninf.chain1$beta,
	gamma_samples_noninf.chain2$beta,
	gamma_samples_noninf.chain3$beta,
	gamma_samples_noninf.chain4$beta)

show_thinning_options(gamma_samples_noninf.chain1$alpha) # 10 is best
show_thinning_options(gamma_samples_noninf.chain1$beta) # 10 is best

# 4. Grab Gamma parameters

final_gamma <- NULL
final_gamma$alpha <- c(
	thin(gamma_samples_noninf.chain1$alpha, 10),
	thin(gamma_samples_noninf.chain1$alpha, 10),
	thin(gamma_samples_noninf.chain1$alpha, 10),
	thin(gamma_samples_noninf.chain1$alpha, 10))
final_gamma$beta <- c(
	thin(gamma_samples_noninf.chain1$beta, 10),
	thin(gamma_samples_noninf.chain1$beta, 10),
	thin(gamma_samples_noninf.chain1$beta, 10),
	thin(gamma_samples_noninf.chain1$beta, 10))

posterior_summary(final_gamma$alpha)
posterior_summary(final_gamma$beta)

# INFORMATIVE
# 1. Tune with Gamma proposal
# Not necessary

# 2. Sample 
gamma_samples_inf.chain1 <- gammaSamp(
	seed = 1, data = all_storms, B = 20000,
	a1 = 60, b1 = 63,
	p=10, q=10, r=10, s=10,
	alpha_start = 0.7, beta_start = 0.2)

gamma_samples_inf.chain2 <- gammaSamp(
	seed = 2, data = all_storms, B = 20000,
	a1 = 60, b1 = 63,
	p=10, q=10, r=10, s=10,
	alpha_start = 1.1, beta_start = 0.3)

gamma_samples_inf.chain3 <- gammaSamp(
	seed = 3, data = all_storms, B = 20000,
	a1 = 60, b1 = 63,
	p=10, q=10, r=10, s=10,
	alpha_start = 0.9, beta_start = 0.4)

gamma_samples_inf.chain4 <- gammaSamp(
	seed = 4, data = all_storms, B = 20000,
	a1 = 60, b1 = 63,
	p=10, q=10, r=10, s=10,
	alpha_start = 0.8, beta_start = 0.5)

# 3. Check for convergence
mcmcplot2(gamma_samples_inf.chain1$alpha)
mcmcplot2(gamma_samples_inf.chain2$alpha)
mcmcplot2(gamma_samples_inf.chain3$alpha)
mcmcplot2(gamma_samples_inf.chain4$alpha)

gelman_rubin(
	gamma_samples_inf.chain1$alpha,
	gamma_samples_inf.chain2$alpha,
	gamma_samples_inf.chain3$alpha,
	gamma_samples_inf.chain4$alpha)

mcmcplot2(gamma_samples_inf.chain1$beta)
mcmcplot2(gamma_samples_inf.chain2$beta)
mcmcplot2(gamma_samples_inf.chain3$beta)
mcmcplot2(gamma_samples_inf.chain4$beta)

gelman_rubin(
	gamma_samples_inf.chain1$beta,
	gamma_samples_inf.chain2$beta,
	gamma_samples_inf.chain3$beta,
	gamma_samples_inf.chain4$beta)

show_thinning_options(gamma_samples_inf.chain1$alpha) # 10 is best
show_thinning_options(gamma_samples_inf.chain1$beta) # 10 is best

# 4. Grab Gamma parameters
posterior_summary(c(
	thin(gamma_samples_inf.chain1$alpha, 10),
	thin(gamma_samples_inf.chain1$alpha, 10),
	thin(gamma_samples_inf.chain1$alpha, 10),
	thin(gamma_samples_inf.chain1$alpha, 10)))

posterior_summary(c(
	thin(gamma_samples_inf.chain1$beta, 10),
	thin(gamma_samples_inf.chain1$beta, 10),
	thin(gamma_samples_inf.chain1$beta, 10),
	thin(gamma_samples_inf.chain1$beta, 10)))


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
weibull_samples.chain1 <- weibullSamp(
	seed = 1, data = all_storms, B = 20000,
	a1 = 60, b1 = 67,
	theta_start = 0.8, lambda_start = 2)

weibull_samples.chain2 <- weibullSamp(
	seed = 2, data = all_storms, B = 20000,
	a1 = 60, b1 = 67,
	theta_start = 1.1, lambda_start = 4)

weibull_samples.chain3 <- weibullSamp(
	seed = 3, data = all_storms, B = 20000,
	a1 = 60, b1 = 67,
	theta_start = 0.9, lambda_start = 3.5)

weibull_samples.chain4 <- weibullSamp(
	seed = 4, data = all_storms, B = 20000,
	a1 = 60, b1 = 67,
	theta_start = 0.9, lambda_start = 2.5)

# 3. Check for convergence
mcmcplot2(weibull_samples.chain1$theta)
mcmcplot2(weibull_samples.chain2$theta)
mcmcplot2(weibull_samples.chain3$theta)
mcmcplot2(weibull_samples.chain4$theta)

gelman_rubin(
	weibull_samples.chain1$theta,
	weibull_samples.chain2$theta,
	weibull_samples.chain3$theta,
	weibull_samples.chain4$theta)

mcmcplot2(weibull_samples.chain1$lambda)
mcmcplot2(weibull_samples.chain2$lambda)
mcmcplot2(weibull_samples.chain3$lambda)
mcmcplot2(weibull_samples.chain4$lambda)

gelman_rubin(
	weibull_samples.chain1$lambda,
	weibull_samples.chain2$lambda,
	weibull_samples.chain3$lambda,
	weibull_samples.chain4$lambda)

show_thinning_options(weibull_samples.chain1$theta) # 5 is best
show_thinning_options(weibull_samples.chain1$lambda) # 5 is best

# 4. Grab Weibull parameters

final_weibull <- NULL
final_weibull$theta <- c(
	thin(weibull_samples.chain1$theta, 5),
	thin(weibull_samples.chain1$theta, 5),
	thin(weibull_samples.chain1$theta, 5),
	thin(weibull_samples.chain1$theta, 5))
final_weibull$lambda <- c(
	thin(weibull_samples.chain1$lambda, 5),
	thin(weibull_samples.chain1$lambda, 5),
	thin(weibull_samples.chain1$lambda, 5),
	thin(weibull_samples.chain1$lambda, 5))

posterior_summary(final_weibull$theta)
posterior_summary(final_weibull$lambda)

###############
# C. Decide on best likelihood (probably tie between weibull, gamma)
###############

# Graphical comparison to empirical density
plot(density(all_storms, from=0), lwd=2, ylim=c(0, 0.4), xlim=c(0,20),
		 main="Model Selection", xlab="Days Elapsed")
legend(10, 0.4, 
			 legend=c("Empirical", "Log-normal", "Gamma", "Weibull"),
		 	 col=c("black", "green", "blue", "red"), lty=c(2,2,2,2), lwd=c(1,2,2,2), cex=1)
curve(dlnorm(x, 0.434, sqrt(1.878)), col="green", lwd=2, lty=2, add=TRUE)
curve(dgamma(x, 0.901, 0.301), col="blue", lwd=2, lty=2, add=TRUE)
curve(dweibull(x, 0.908, 2.877), col="red", lwd=2, lty=2, add=TRUE)

# Show replicated moment distribution with original data
show_replicate_analysis(all_storms, rlnorm, final_lognorm$mu, sqrt(final_lognorm$sig2)) # remember sqrt!
show_replicate_analysis(all_storms, rgamma, final_gamma$alpha, final_gamma$beta)
show_replicate_analysis(all_storms, rweibull, final_weibull$theta, final_weibull$lambda)
dev.off()

# Decide based on DIC
calc_DIC(
	data = all_storms, 
	ddistribution = dlnorm,
	params = cbind(final_lognorm$mu, sqrt(final_lognorm$sig2)))

calc_DIC(
	data = all_storms, 
	ddistribution = dgamma,
	params = cbind(final_gamma$alpha, final_gamma$beta))

calc_DIC(
	data = all_storms, 
	ddistribution = dweibull,
	params = cbind(final_weibull$theta, final_weibull$lambda))

###############
# D. Compare credible intervals of pre- and post-election 
# (there is a difference, but not stat significant because intervals overlap)
###############

# Use Weibull because paramters have interesting interpretation
source("weibull_sampler.R")

# Train using pre-election data
pre_weibull_samples.chain1 <- weibullSamp(
	seed = 1, data = pre_election, B = 20000,
	a1 = 60, b1 = 67,
	theta_start = 0.8, lambda_start = 2)
pre_weibull_samples.chain2 <- weibullSamp(
	seed = 2, data = pre_election, B = 20000,
	a1 = 60, b1 = 67,
	theta_start = 1.1, lambda_start = 4)
pre_weibull_samples.chain3 <- weibullSamp(
	seed = 3, data = pre_election, B = 20000,
	a1 = 60, b1 = 67,
	theta_start = 0.9, lambda_start = 3.5)
pre_weibull_samples.chain4 <- weibullSamp(
	seed = 4, data = pre_election, B = 20000,
	a1 = 60, b1 = 67,
	theta_start = 0.9, lambda_start = 2.5)

mcmcplot2(pre_weibull_samples.chain1$theta)
mcmcplot2(pre_weibull_samples.chain2$theta)
mcmcplot2(pre_weibull_samples.chain3$theta)
mcmcplot2(pre_weibull_samples.chain4$theta)

gelman_rubin(
	pre_weibull_samples.chain1$theta,
	pre_weibull_samples.chain2$theta,
	pre_weibull_samples.chain3$theta,
	pre_weibull_samples.chain4$theta)

mcmcplot2(pre_weibull_samples.chain1$lambda)
mcmcplot2(pre_weibull_samples.chain2$lambda)
mcmcplot2(pre_weibull_samples.chain3$lambda)
mcmcplot2(pre_weibull_samples.chain4$lambda)

gelman_rubin(
	pre_weibull_samples.chain1$lambda,
	pre_weibull_samples.chain2$lambda,
	pre_weibull_samples.chain3$lambda,
	pre_weibull_samples.chain4$lambda)

show_thinning_options(pre_weibull_samples.chain1$theta) # 5 is best
show_thinning_options(pre_weibull_samples.chain1$lambda) # 5 is best

pre_weibull <- NULL
pre_weibull$theta <- c(
	thin(pre_weibull_samples.chain1$theta, 5),
	thin(pre_weibull_samples.chain1$theta, 5),
	thin(pre_weibull_samples.chain1$theta, 5),
	thin(pre_weibull_samples.chain1$theta, 5))
pre_weibull$lambda <- c(
	thin(pre_weibull_samples.chain1$lambda, 5),
	thin(pre_weibull_samples.chain1$lambda, 5),
	thin(pre_weibull_samples.chain1$lambda, 5),
	thin(pre_weibull_samples.chain1$lambda, 5))

posterior_summary(pre_weibull$theta)
posterior_summary(pre_weibull$lambda)

# Train using post-election data
post_weibull_samples.chain1 <- weibullSamp(
	seed = 1, data = post_election, B = 20000,
	a1 = 60, b1 = 67,
	theta_start = 0.8, lambda_start = 2)
post_weibull_samples.chain2 <- weibullSamp(
	seed = 2, data = post_election, B = 20000,
	a1 = 60, b1 = 67,
	theta_start = 1.1, lambda_start = 4)
post_weibull_samples.chain3 <- weibullSamp(
	seed = 3, data = post_election, B = 20000,
	a1 = 60, b1 = 67,
	theta_start = 0.9, lambda_start = 3.5)
post_weibull_samples.chain4 <- weibullSamp(
	seed = 4, data = post_election, B = 20000,
	a1 = 60, b1 = 67,
	theta_start = 0.9, lambda_start = 2.5)

mcmcplot2(post_weibull_samples.chain1$theta)
mcmcplot2(post_weibull_samples.chain2$theta)
mcmcplot2(post_weibull_samples.chain3$theta)
mcmcplot2(post_weibull_samples.chain4$theta)

gelman_rubin(
	post_weibull_samples.chain1$theta,
	post_weibull_samples.chain2$theta,
	post_weibull_samples.chain3$theta,
	post_weibull_samples.chain4$theta)

mcmcplot2(post_weibull_samples.chain1$lambda)
mcmcplot2(post_weibull_samples.chain2$lambda)
mcmcplot2(post_weibull_samples.chain3$lambda)
mcmcplot2(post_weibull_samples.chain4$lambda)

gelman_rubin(
	post_weibull_samples.chain1$lambda,
	post_weibull_samples.chain2$lambda,
	post_weibull_samples.chain3$lambda,
	post_weibull_samples.chain4$lambda)

show_thinning_options(post_weibull_samples.chain1$theta) # 5 is best
show_thinning_options(post_weibull_samples.chain1$lambda) # 5 is best

post_weibull <- NULL
post_weibull$theta <- c(
	thin(post_weibull_samples.chain1$theta, 5),
	thin(post_weibull_samples.chain1$theta, 5),
	thin(post_weibull_samples.chain1$theta, 5),
	thin(post_weibull_samples.chain1$theta, 5))
post_weibull$lambda <- c(
	thin(post_weibull_samples.chain1$lambda, 5),
	thin(post_weibull_samples.chain1$lambda, 5),
	thin(post_weibull_samples.chain1$lambda, 5),
	thin(post_weibull_samples.chain1$lambda, 5))

posterior_summary(post_weibull$theta)
posterior_summary(post_weibull$lambda)

# Look at median of distribution via parameters
pre_median <- pre_weibull$lambda * log(2)^(1/pre_weibull$theta)
post_median <- post_weibull$lambda * log(2)^(1/post_weibull$theta)

posterior_summary(pre_median)
posterior_summary(post_median)

caterplot2(pre_median, "Pre-\nElection", post_median, "Post-\nElection")
title(main="Median Days Elapsed\nSince Last Tweetstorm", 
			xlab="Time (days)")

# Look at Prob(> 1 week silence)
n_pre <- post_weibull$theta %>% length
n_post <- post_weibull$theta %>% length
mean(rweibull(n_pre, pre_weibull$theta, pre_weibull$lambda) > 7)
mean(rweibull(n_post, post_weibull$theta, post_weibull$lambda) > 7)

###############
# E. How well does pre-election model/data do for post-election realized data?
###############

plot(density(post_election, from=0), lwd=2, col="black", xlim=c(0, 20),
		 main="Pre-election model vs.\npost-election data")
curve(dweibull(x, 0.8359, 3.4118), add=TRUE, col="purple", lty=2, lwd=2)
legend(6, 0.2, 
			 legend=c("Empirical (post-election)", "Model (pre-election)"),
			 col=c("black", "purple", "blue", "red"), lty=c(1,2), lwd=c(2,2), cex=1)


