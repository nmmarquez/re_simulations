rm(list=ls())
library(MASS)
set.seed(1)

simulate_time_series <- function(years, draws, corr){
    X <- matrix(runif(draws*years), nrow=years, ncol=draws)
    corr_mat <- corr**abs(outer(0:(years-1), 0:(years-1), "-"))
    mvdat <- t(mvrnorm(n=draws, mu=0 * 1:years, Sigma=corr_mat, empirical=TRUE))
    ranks <- t(apply(mvdat, 1, rank, ties.method="first"))
    sorted_X <- t(apply(X, 1, sort))
    t(sapply(1:years, function(x) sorted_X[x,][ranks[x,]]))
}

draws <- 1000
years <- 20
corr <- .98

corr_X <- simulate_time_series(years, draws, corr)

sapply(2:years, function(x) cor(corr_X[x,], corr_X[x-1,]))

ages <- 10
age_corr <- .85

age_simulations <- lapply(1:ages, function(x) 
    simulate_time_series(years, draws, corr))

X_age <- t(sapply(age_simulations, colSums))
corr_mat_age <- age_corr**abs(outer(0:(ages-1), 0:(ages-1), "-"))
mvdat <- t(mvrnorm(n=draws, mu=0 * 1:ages, Sigma=corr_mat_age, empirical=TRUE))
ranks <- t(apply(mvdat, 1, rank, ties.method="first"))
rank_data_sum <- t(apply(X_age, 1, rank, ties.method="first"))
sorted_X_age <- t(apply(X_age, 1, sort))
sorted_age_sims <- lapply(1:ages, function(x) age_simulations[[x]])

lapply(1:ages, 1, function(x) age_simulations[[x]][,])

t(sapply(1:years, function(x) sorted_X[x,][ranks[x,]]))