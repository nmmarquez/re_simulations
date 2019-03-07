rm(list=ls())
library(ar.matrix)
library(ggplot2)
library(dplyr)
library(tidyr)

M <- 50
N <- 1000
sigma.ar <- 1
rho.ar <- .9
margSD <- sqrt(sigma.ar^2 / (1 - rho.ar^2))

Q <- Q.AR1(M, sigma.ar, rho.ar)
Sigma <- (sigma.ar^2/(1-rho.ar^2)) * rho.ar^abs(
    matrix(rep(1:M, M), nrow=M) - matrix(rep(1:M, each=M), nrow=M))

sim1 <- r.AR1(N, M, sigma.ar, rho.ar) 
apply(sim1, 1, sum)

sim2 <- t(apply(sim1, 1, function(x){
    ones <- matrix(rep(1, M), 1, M)
    x - Sigma %*% t(ones) %*% (ones %*% Sigma %*% t(ones))^-1 %*% (ones %*% x)
}))
apply(sim2, 1, sum)

apply(sim1, 2, sd)
apply(sim2, 2, sd)

# try with just regular iid 
M <- 50
N <- 1000
Q <- diag(seq(1, 5, length.out=M)^-1)
Sigma <- diag(seq(1, 5, length.out=M))
sim1 <- sim.AR(N, Q)
head(apply(sim1, 1, sum))

sim2 <- t(apply(sim1, 1, function(x){
    ones <- matrix(rep(1, M), 1, M)
    x - Sigma %*% t(ones) %*% (ones %*% Sigma %*% t(ones))^-1 %*% (ones %*% x)
}))
head(apply(sim2, 1, sum))

data.frame(
    sim1=apply(sim1, 2, sd),
    sim2=apply(sim2, 2, sd),
    true=seq(1, 5, length.out=M)^.5)

# something seems off so lets try a new experiment where we check bias of estimator
M <- 50
N <- 10000
sigma.iid <- sqrt(100)
Q <- diag(rep(sigma.iid^-2, M))
Sigma <- diag(rep(sigma.iid^-2, M))
sim1 <- sim.AR(N, Q)
head(apply(sim1, 1, sum))

sim2 <- t(apply(sim1, 1, function(x){
    ones <- matrix(rep(1, M), 1, M)
    x - Sigma %*% t(ones) %*% (ones %*% Sigma %*% t(ones))^-1 %*% (ones %*% x)
}))
head(apply(sim2, 1, sum))

data.frame(sim=apply(sim1, 2, sd), simstar=apply(sim2, 2, sd)) %>%
    gather(key="Simulation", value="SD") %>%
    ggplot(aes(x=SD, fill=Simulation, group=Simulation)) +
    geom_density(alpha=.2) +
    geom_vline(xintercept=sigma.iid, linetype=3) +
    theme_classic()

