rm(list=ls())
library(tidyr)
library(dplyr)
library(Matrix)
library(mvtnorm)
library(sparseMVN)
library(ggplot2)

set.seed(123)
# parameters of simulation
M <- 100
N <- 5000
sigma.ar <- 1
rho.ar <- .9

# define a Q(precision matrix) for AR1
Q.AR1 <- function (M, sigma, rho, sparse = FALSE, vcov = FALSE){
    if (sigma <= 0) 
        stop("sigma paramter must be greater than 0.")
    Q <- matrix(0, nrow = M, ncol = M)
    Q[1, 1] <- 1
    for (i in 2:M) {
        Q[i, i] <- 1 + rho^2
        Q[i - 1, i] <- -1 * rho
        Q[i, i - 1] <- -1 * rho
    }
    Q[M, M] <- 1
    Q <- (1/sigma^2) * Q
    if (vcov) 
        Q <- solve(Q)
    if (sparse) 
        Q <- Matrix(Q, sparse = TRUE)
    Q
}

# cholesky simulation
sim.AR <- function (n, Q){
    library(Matrix)
    library(sparseMVN)
    rmvn.sparse(n, rep(0, nrow(Q)), Cholesky(Matrix(Q, sparse = TRUE)), T)
}

# calculate marginal sd
margSD <- sqrt(sigma.ar^2 / (1 - rho.ar^2))

# simulate by hand

handAR <- t(sapply(1:N, function(i){
    vec <- rep(0, M)
    for(j in 2:M){
        vec[j] <- rnorm(1, rho.ar*vec[j-1], sigma.ar)
    }
    vec
}))

# arima.sim
statsAR1 <- t(sapply(1:N, function(i){
    vec <- arima.sim(n=M, list(ar=rho.ar))
    vec
}))

# cholesky decomp of Q manually
cholAR <- sim.AR(N, Q.AR1(M, sigma.ar, rho.ar, T))
cholARReg <- cholAR - cholAR[,1]

mvnAR <- rmvnorm(N, rep(0,M), (sigma.ar^2/(1-rho.ar^2)) * rho.ar^abs(
    matrix(rep(1:M, M), nrow=M) - matrix(rep(1:M, each=M), nrow=M)))
mvnARReg <- mvnAR - mvnAR[,1]

data.frame(
    hand=apply(handAR, 2, sd),
    stats=apply(statsAR1, 2, sd),
    cholQ=apply(cholAR, 2, sd),
    mvn=apply(mvnAR, 2, sd),
    Time=1:M) %>%
    gather(key="Model", value="SD", -Time) %>%
    ggplot(aes(Time, SD, group=Model, color=Model)) +
    geom_line() +
    theme_classic() +
    geom_hline(yintercept=margSD, linetype=2, color="red")
