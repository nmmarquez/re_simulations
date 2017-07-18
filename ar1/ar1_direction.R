rm(list=ls())
set.seed(123)
pacman::p_load(clusterPower, mvtnorm, LaplacesDemon)

# simulate an ar1 with a for loop
sim_ar1 <- function(N, rho, sigma){
    timeseries <- rnorm(N, 0, sd=sigma)
    for(i in 2:N){
        timeseries[i] <- rho*timeseries[i-1] + timeseries[i]
    }
    timeseries
}

# set the parameters
N <- 1000
rho <- .95
sigma <- 1.
ts_ <- sim_ar1(N, rho, sigma)
plot(ts_, type="l")

# utility function for trnasformations
transforms <- function(x){
    c(rho=expit(x[1]), sigma=exp(x[2]))
}

# evaluate the nll of an ar1 process
ar1nll1 <- function(params, ts_){
    paramst <- transforms(params)
    rho <- paramst["rho"]
    sigma <- paramst["sigma"]
    N <- length(ts_)
    tshat <- (ts_ * rho)[1:(N-1)]
    tsdat <- ts_[2:N]
    return(sum(-1 * dnorm(tsdat, tshat, sigma, log=TRUE)))
}

# evaluate the nll of an ar1 process using mvtnorm
ar1nll2 <- function(params, ts_){
    paramst <- transforms(params)
    rho <- paramst["rho"]
    sigma <- paramst["sigma"]
    N <- length(ts_)
    Q <- matrix(0, nrow=N, ncol=N)
    Q[1, 1] <- 1. / sigma**2
    for(i in 2:N){
        Q[i, i] <- (1. + rho**2) / sigma**2
        Q[i-1, i] <- (-1. * rho) / sigma**2
        Q[i, i-1] <- (-1. * rho) / sigma**2
    }
    Q[N, N] <- 1. / sigma**2
    return(-1 * dmvnp(ts_, rep(0, N), Q, log=TRUE))
}

# evaluate time series in the original direction
transforms(optim(c(0, 0), ar1nll1, ts_=ts_)$par)
# rho     sigma 
# 0.9427267 0.9915066 

# use the rev command to reverse it then evaluate it
transforms(optim(c(0, 0), ar1nll1, ts_=rev(ts_))$par)
# rho     sigma 
# 0.9423960 0.9912705 

# evaluate time series in the original direction
transforms(optim(c(0, 0), ar1nll2, ts_=ts_)$par)
# rho     sigma 
# 0.9427267 0.9915066 

# use the rev command to reverse it then evaluate it
transforms(optim(c(0, 0), ar1nll2, ts_=rev(ts_))$par)
# rho     sigma 
# 0.9423960 0.9912705 