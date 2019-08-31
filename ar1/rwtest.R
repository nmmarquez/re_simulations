rm(list=ls())
library(tidyverse)

simulateRW1 <- function(N, sigma_=1, burnin=1000){
    resultsplus <- cumsum(rnorm(burnin+N, 0, sigma_))
    resultsT <- resultsplus[-1:-burnin]
    results <- resultsT - mean(resultsT)
    return(results)
}

simulateRW2 <- function(N, sigma_=1, burnin=1000){
    # resultsplus <- rnorm(burnin+N, 0, (sigma_/2))
    # for(i in 3:(burnin+N)){
    #     resultsplus[i] <- resultsplus[i] +
    #         ((2/3)*resultsplus[i-1]) + ((1/3)*resultsplus[i-2])
    # }
    resultsplus <- rnorm(burnin+N, 0, (sigma_/2))
    for(i in 3:(burnin+N)){
        resultsplus[i] <- resultsplus[i] +
            ((.5)*resultsplus[i-1]) + ((.5)*resultsplus[i-2])
    }
    resultsT <- resultsplus[-1:-burnin]
    results <- resultsT - mean(resultsT)
    return(results)
}


rw2Q <- function(M, sparse=T){
    Q <- Matrix::Matrix(0, nrow = M, ncol = M)
    Q[1, 1:3] <- c(1, -2, 1)
    Q[2, 1:4] <- c(-2, 5, -4, 1)
    for (i in 3:(M-2)){
        Q[i, (i-2):(i+2)] <- c(1, -4, 6, -4, 1)
    }
    Q[M-1, (M-3):M] <- c(1, -4, 5, -2)
    Q[M, (M-2):M] <- c(1, -2, 1)
    if(sparse)
        Q <- Matrix::Matrix(Q, sparse = TRUE)
    Q
}

testDifferences <- function(x, second=T){
    if(second){
        xdiff <- sapply(3:length(x), function(i){
            x[i] - (2 * x[i-1]) + x[i-2]
        })
    }
    else{
        xdiff <- sapply(2:length(x), function(i){
            x[i] - x[i-1]
        })
    }
    xdiff
}


y <- simulateRW1(1000, 1)
plot(y)
plot(density(testDifferences(y, F)))
quantile(testDifferences(y, F), probs = c(.025, .975))

x <- simulateRW2(400, 2)
x <- as.vector(arima.sim(list(ar=c(.499, .499)), n=400, function(z) rnorm(z, sd=.5)))
plot(x)
plot(density(testDifferences(x)))
quantile(testDifferences(x), probs = c(.025, .975))

# x <- as.vector(arima.sim(list(ar=c(.66, .33)), n=100))
# plot(x)
# plot(density(testDifferences(x)))
# quantile(testDifferences(x), probs = c(.025, .975))

densRW2 <- function(x, sigma_=1, log=T){
    kappa_ <- 1/(sigma_^2)
    N <- length(x)
    Q <- rw2Q(N) * kappa_
    ((N-2)/2) * log(kappa_) + (-.5 * ((t(x) %*% Q %*% x)[1,1]))
}

testDF <- tibble(test = seq(.75, 2.5, by=.05)) %>%
    mutate(est=sapply(test, function(s) -densRW2(x, s)))

testDF %>%
    ggplot(aes(x=test, y=est)) +
    geom_point()

testDF %>%
    arrange(est)
