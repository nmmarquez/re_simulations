# standard deviation vs standard error
rm(list=ls())
set.seed(123)

# true population mean
mu <- 3
# true population standard deviation
sigma <- .2
# sample size
N <- 1000

# given the true parameters of the population take a sample
# we will assume a normal distribution
sample.data <- rnorm(N, mean=mu, sd=sigma)

# this is the estimated population mean from the sample which is a sample statistic
mean(sample.data)

# this is the estimated population sd from the sample
sd(sample.data)

# this is an estimate of the standard error of the mean 
sd(sample.data) / sqrt(N)

# this is an estimand, or true, standard error of the mean 
# which is the standard deviation of the sample mean estimate
sigma / sqrt(N)
