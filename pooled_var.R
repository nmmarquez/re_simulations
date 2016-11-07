rm(list=ls())
set.seed(123)
library(ggplot2)

run_simulation <- function(rdist, N_max=100, ...){
    sample_sizes <- sample(4:N_max, 10, replace=T)
    samples_collected <- sapply(sample_sizes, rdist, ...)
    sample_means <- sapply(samples_collected, mean)
    sample_variances <- sapply(samples_collected, var)
    
    sig_mean <- mean(sample_variances)
    
    # pooling variances
    sig_pool <- (sum((sample_sizes - 1) * sample_variances) / 
                     sum(sample_sizes -1))
    c(sig_mean, sig_pool)
}


run_full_sim <- function(M, xint, rdist, N_max=100, ...){
    estimates <- data.frame(est=c(sapply(1:M, function(x) 
        run_simulation(rdist, N_max, ...))))
    estimates$type <- c(rep("mean_sig", M), rep("pooled_sig", M))
    
    ggplot(estimates, aes(est, fill = type, colour = type)) +
        geom_density(alpha = 0.1) + geom_vline(xintercept=xint) +
        xlab("Estimated Variance")
}

# sample size for alll simulations
M <- 10000

# normal distribution simulation
mu <- 3.2
sigma <- 1.3

# simulation
run_full_sim(M, sigma**2, rnorm, mean=mu, sd=sigma)

# binomial parameters
p <- .3
run_full_sim(M, p*(1-p), rbinom, size=1, prob=p)

# small binomial parameters
p_small <- .01
run_full_sim(M, p_small*(1-p_small), rbinom, size=1, prob=p_small)

# poisson sim
lambda <- 2
run_full_sim(M, lambda, rpois, lambda=lambda)

# small
lambda_small <- .2
run_full_sim(M, lambda_small, rpois, lambda=lambda_small)
