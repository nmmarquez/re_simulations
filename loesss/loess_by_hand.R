rm(list=ls())
pacman::p_load(ggplot2, data.table)
set.seed(123)

# simulate some data with a sin like relationship
x <- runif(500, min=-6, max=6)
y <- sin(x) + rnorm(length(x))

DT <- data.table(x=x, y=y)
ggplot(data=DT, aes(x=x, y=y)) + geom_point() + geom_smooth(method=lm, se=FALSE)


loess_ <- function(y, x, alpha=.75, degree=1){
    N <- length(x)
    dist <- sapply(x, function(z) abs(x - z))
    # find the alpha percentile highest distance
    cutoffarray <- apply(dist, 1, quantile, probs=alpha)
    # calculate weights base on origin lab page above
    cuberootweights <- 1 - (dist / cutoffarray)**3
    cuberootweights[cuberootweights < 0] <- 0
    weights <- cuberootweights**3
    # run a WLS for each point based on wiki page above
    X <- sapply(0:degree, function(i) x**i)
    Y <-  matrix(y, length(y), 1)
    betas <- X * 0
    for(i in 1:N){
        W <- diag(weights[, i])
        betas[i,] <- t(solve(t(X) %*% W %*% X) %*% (t(X) %*% W %*% Y))
    }
    preds <- rowSums(X * betas)
    return(preds)
}

DThat <- data.table(x=x, y=loess_(y, x, alpha=.5, degree=2))

ggplot(data=DT, aes(x=x, y=y)) + geom_point() + 
  geom_line(data=DThat, color="red", size=2)
