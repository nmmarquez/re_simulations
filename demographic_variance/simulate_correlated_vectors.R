rm(list=ls())
library(MASS)
set.seed(1)

simulate_time_series <- function(years, draws, corr, print=FALSE){
    X <- matrix(runif(draws*years), nrow=years, ncol=draws)
    if(print){
        print(sapply(2:years, function(x) cor(X[x,], X[x-1,])))
    }
    corr_mat <- corr**abs(outer(0:(years-1), 0:(years-1), "-"))
    mvdat <- t(mvrnorm(n=draws, mu=0 * 1:years, Sigma=corr_mat, empirical=TRUE))
    ranks <- t(apply(mvdat, 1, rank, ties.method="first"))
    sorted_X <- t(apply(X, 1, sort))
    t(sapply(1:years, function(x) sorted_X[x,][ranks[x,]]))
}

draws <- 10
years <- 5
corr <- .98

corr_X <- simulate_time_series(years, draws, corr, print=TRUE)
sapply(2:years, function(x) cor(corr_X[x,], corr_X[x-1,]))

ages <- 20
age_corr <- .85

age_simulations <- lapply(1:ages, function(x) 
    simulate_time_series(years, draws, corr))

sapply(1:ages, function(x) sapply(2:years, function(y)
    cor(age_simulations[[x]][y,], age_simulations[[x]][y-1,])))

sapply(2:ages, function(x) sapply(1:years, function(y)
    cor(age_simulations[[x]][y,], age_simulations[[x-1]][y,])))

X_age <- t(sapply(age_simulations, colSums))
corr_mat_age <- age_corr**abs(outer(0:(ages-1), 0:(ages-1), "-"))
mvdat <- t(mvrnorm(n=draws, mu=0 * 1:ages, Sigma=corr_mat_age, empirical=TRUE))
ranks <- t(apply(mvdat, 1, rank, ties.method="first"))
sorted_X_age <- t(apply(X_age, 1, function(x) sort(x, index.return=TRUE)$ix))
sorted_age_sims <- lapply(1:ages, function(x) 
    age_simulations[[x]][,sorted_X_age[x,]])

age_time_corr_X <- lapply(1:ages, function(x)
    sorted_age_sims[[x]][,ranks[x,]])

sapply(1:ages, function(x) sapply(2:years, function(y)
    cor(age_time_corr_X[[x]][y,], age_time_corr_X[[x]][y-1,])))

sapply(2:ages, function(x) sapply(1:years, function(y)
    cor(age_time_corr_X[[x]][y,], age_time_corr_X[[x-1]][y,])))

