rm(list=ls())
library(ar.matrix)

M <- 320
N <- 1000
sigma.ar <- 2
rho.ar <- .9
Sigma <- Q.iid(M, sigma.ar^-1)

X <- r.iid(N, M, sigma.ar)

summary(apply(X, 1, sum))
summary(apply(X, 1, sd))

Xstar <- t(apply(X, 1, function(x){
    A <- matrix(rep(1, M), nrow=1, ncol=M)
    x - Sigma %*% t(A) %*% (A %*% Sigma %*% t(A))^-1 %*% (A %*% x)
}))

sapply(1:N, function(i){
    cor(X[i,], Xstar[i,])
})

summary(apply(Xstar, 1, sum))
summary(apply(Xstar, 1, sd))
