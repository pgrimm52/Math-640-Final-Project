
prior <- function(alpha, beta, p, q, r, s){
	beta^(alpha*s) / gamma(alpha)^r * p^(alpha-1) * exp(-beta*q)
}

curve(dgamma(x, 7, 2), xlim=c(0, 10))


x1 <- seq(0.1, 1, by=.1)
x2 <- x1
dens <- matrix(prior(
	alpha = expand.grid(x1, x2)[, 1],
	beta = expand.grid(x1, x2)[, 2],
	p=10, q=10, r=10, s=10), ncol=length(x1))

plotly::plot_ly(z = ~ dens) %>% add_surface()

###

small_data <- all_storms[sample(1:length(all_storms), 20)]
small_data <- all_storms

gamma_samples <- gammaSamp(data = small_data, B = 10000,
													 a1 = 64, b1 = 67,
													 p=10, q=10, r=10, s=10,
													 alpha_start = 1, beta_start = 1)
lapply(gamma_samples, mean)

hist(small_data)
curve(10*dgamma(x, 1.3, 0.9), add=TRUE)

gamma_samples <- gammaSamp(data = small_data, B = 10000,
													 a1 = 64, b1 = 67,
													 p=1, q=0, r=0, s=0,
													 alpha_start = 1, beta_start = 1)
lapply(gamma_samples, mean)
