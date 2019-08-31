rm(list=ls())
set.seed(123)
library(tidyverse)
library(TMB)
library(ar.matrix)
library(mvtnorm)
library(sparseMVN)
library(Matrix)

simulateRW2 <- function(N, sigma_=1, burnin=1000){
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

trueValue <- 1.4

x <- simulateRW2(400, trueValue)
plot(x)
plot(density(testDifferences(x)))
quantile(testDifferences(x), probs = c(.025, .975))

densRW2 <- function(log_sigma, x_=x, log=T, alt=F, offset=.00001){
    sigma_ <- exp(log_sigma)
    kappa_ <- 1/(sigma_^2)
    N <- length(x_)
    Q <- rw2Q(N) * kappa_
    nll <- -1 * (((N-2)/2) * log(kappa_) + (-.5 * ((t(x_) %*% Q %*% x_)[1,1])))
    if(!log){
        nll <- exp(nll)
    }
    if(alt){
        diag(Q) <- diag(Q) + .00001
        CH <- Cholesky(Q)
        nll <- -dmvn.sparse(x_, rep(0, nrow(Q)), CH)
    }
    return(nll)
}

modelRun <- function(x, verbose=T){
    model <- "rw2"
    compile(paste0(model, ".cpp"))
    dyn.load(dynlib(model))
    Data <- list(x=x)
    Params <- list(log_sigma=0)
    
    Obj <- MakeADFun(data=Data, parameters=Params, DLL=model, silent=!verbose)
    
    Opt <- nlminb(start=Obj$par, objective=Obj$fn, 
                  gradient=Obj$gr,
                  control=list(eval.max=1e6, iter.max=1e6))
    
    sdrep <- sdreport(Obj, getJointPrecision=TRUE)
    dyn.unload(dynlib(model))
    
    list(obj=Obj, opt=Opt, sd=sdrep)
}


test <- modelRun(x)
test2 <- optim(
    0, densRW2, method="Brent", lower=-100, upper=100, hessian=T, alt=T)


data.frame(
    method=c("tmb", "optim"),
    mu=c(unname(exp(test$sd$par.fixed)), exp(test2$par)),
    lwr=exp(c(
        test$sd$par.fixed - 1.96 * test$sd$cov.fixed^.5,
        test2$par - 1.96 * test2$hessian^-.5)),
    upr=exp(c(
        test$sd$par.fixed + 1.96 * test$sd$cov.fixed^.5,
        test2$par + 1.96 * test2$hessian^-.5)),
    true=c(trueValue, trueValue))
