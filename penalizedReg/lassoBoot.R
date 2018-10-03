rm(list=ls())
library(MASS)
library(glmnet)
library(dplyr)
library(ggplot2)

N <- 90 # number of observations
betas <- c(-5, 4, -3, 2, -1, rep(0, 95)) # betas with lots of zeros
X <- cbind(rep(1, N), matrix(rnorm(N*(length(betas)-1)), nrow=N)) # covariates
sigma.epsilon <- .4 # random error
Y <- rnorm(N, c(X %*% betas), sigma.epsilon) # observed values

# lasso bootstrap
lassoBetaBoot <- function(X, Y, boots=100, lambdas=10^seq(5, -5, length=50), ...){
    if(nrow(X) != length(Y)){
        stop("Number of covariate observations does not match outcome.")
    }
    N <- nrow(X) # number of obs
    sapply(1:boots, function(i){
        bsamp <- sample.int(N, replace=T) # sample with replacement
        Xboot <- X[bsamp,] 
        Yboot <- Y[bsamp]
        # fit full data
        fitboot <- glmnet(Xboot, Yboot, alpha=1, lambda=lambdas)
        # fit with cross validation to find best lamda
        lambmin <- cv.glmnet(Xboot, Yboot, alpha=1, lambda=lambdas)$lambda.min
        # predict coefficients with best lambda
        as.vector(predict(fitboot, s=lambmin, exact=T, type = 'coefficients'))
    })
    
}

boots <- 100
betaResults <- lassoBetaBoot(X[,-1], Y, boots)
data.frame(
    beta=rep(1:length(betas), boots),
    draw=rep(1:boots, each=length(betas)),
    est=c(betaResults)) %>%
    group_by(beta) %>%
    summarize(
        mu=median(est), 
        lwr=quantile(est, probs=.025),
        upr=quantile(est, probs=.975)) %>%
    mutate(`True Effect`=beta<=5) %>%
    ggplot(aes(beta, mu, ymin=lwr, ymax=upr, color=`True Effect`)) +
    geom_errorbar() +
    theme_classic()
