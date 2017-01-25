rm(list=ls())
library(MASS)

N <- 1000
obs <- 3
data <- sapply(1:obs, function(x) rnorm(N))
corr <- .7
cor(data)

copulate_data <- function(X, corr){
    if(length(corr) == 1){
        corr <- matrix(data=corr, nrow=ncol(X), ncol=ncol(X))
        diag(corr) <- 1
    }
    mvdat <- mvrnorm(n=nrow(X), mu=0 * 1:ncol(X), Sigma=corr, empirical=T)
    ranks <- apply(mvdat, 2, rank, ties.method="first")
    sorted_X <- apply(X, 2, sort)
    sapply(1:ncol(X), function(x) sorted_X[,x][ranks[,x]])
}

new_data <- copulate_data(data, corr)
cor(new_data)
