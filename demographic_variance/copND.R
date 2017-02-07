rm(list=ls())
library(MASS)
set.seed(1)

# Dunction for creating AR draws each draw has an AR process for x years
simulate_time_series <- function(years, draws, corr, print=FALSE){
    X <- matrix(runif(draws*years), nrow=years, ncol=draws)
    corr_mat <- corr**abs(outer(0:(years-1), 0:(years-1), "-"))
    mvdat <- t(mvrnorm(n=draws, mu=0 * 1:years, Sigma=corr_mat, empirical=TRUE))
    ranks <- t(apply(mvdat, 1, rank, ties.method="first"))
    sorted_X <- t(apply(X, 1, sort))
    t(sapply(1:years, function(x) sorted_X[x,][ranks[x,]]))
}

draws <- 100 # number of draws
years <- 20 # number of years
corr <- .98 # correlation over time
ages <- 10 # number of age groups we have
age_corr <- .85 # correlation we want between adjacent age groups

# simulate draws of AR time series for independent age groups
age_simulations <- lapply(1:ages, function(x) 
    simulate_time_series(years, draws, corr))

# restructure data so it is in 3D array with dims being c(time, age, draws)
Xtad <- aperm(array(c(sapply(age_simulations, function(x) x)), 
                    dim=c(years, draws,ages)), c(1, 3, 2))

# check out the correlation over time it looks good thanks to first function
sapply(1:ages, function(x) sapply(2:years, function(y)
    cor(Xtad[y,x,], Xtad[y-1,x,])))

# since ages were created independently however they are uncorrelated
sapply(2:ages, function(x) sapply(1:years, function(y)
    cor(Xtad[y,x,], Xtad[y,x-1,])))

# here I am going to devise the corr mat I want for ages but this is most 
# likely going to be derived from data for other cases
corr_mat <- age_corr**abs(outer(0:(ages-1), 0:(ages-1), "-"))

# this is the function that takes a 3D array and leaves the first dimension (in
# our test case that dimension is time) unchnaged while sorting the 
# 3rd dimension (this is probably always gonna be the draws dimension) 
# in order to get the desired correlation in the 2nd dimension (for us it is
# age but it could just as easily be country or cause or SDI component whatever)
draw2Dcopula <- function(X, cor_mat){
    L <- dim(X)[2]
    D <- dim(X)[3]
    Xsum <- apply(X, c(2, 3), sum)
    mvdat <- mvrnorm(n=D, mu=0 * 1:L, Sigma=cor_mat, empirical=TRUE)
    ranks <- apply(mvdat, 2, rank, ties.method="first")
    sortedXsim <- apply(Xsum, 1, function(x) sort(x, index.return=TRUE)$ix)
    sortedX <- Xtad
    for(i in 1:L){
        sortedX[,i,] <- Xtad[,i,sortedXsim[,i]]
    }
    Xcorr <- sortedX
    for(i in 1:L){
        Xcorr[,i,] <- sortedX[,i,ranks[,i]]
    }
    Xcorr
}

# lets 2D copulate the data now
Xcorr <- draw2Dcopula(Xtad, corr_mat)

# correlation over our first dimension (time) remains unchanged
sapply(1:ages, function(x) sapply(2:years, function(y)
    cor(Xcorr[y,x,], Xcorr[y-1,x,])))

# correlation over our second dimension (age) now approaches desired target
sapply(2:ages, function(x) sapply(1:years, function(y)
    cor(Xcorr[y,x,], Xcorr[y,x-1,])))
