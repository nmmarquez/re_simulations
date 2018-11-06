set.seed(123)
# library to run multinomial regression
library(nnet)
library(MASS)
library(mvtnorm)
rm(list=ls())

N <- 10000 # number of observations
M <- 2 # number of non-intercept covariates
G <- 5 # number of groups

# generate the covariates
covs <- sapply(1:M, function(x) rnorm(N))
covs <- cbind(rep(1, N), covs)

# generate the betas
betas <- sapply(1:(G-1), function(x) runif(M+1, min = -.7, max=.4))

# the response var in log_ratio_space
ratio_ref <- exp(covs %*% betas)

# the response var in probability space
p1 <- apply(ratio_ref, 1, function(x) 1 / (1 + sum(x)))
y_prob <- cbind(p1, ratio_ref * p1)

# we will add some noise to the process which will be the noise associted
# wuith being able to calculate the true underlying probability association
# process. Our response variable will then be a matrix of N observations
# with a number of columns equal to the number of groups and a sum to one
# row-wise constraint.
clean <- 10 # how "clean" is the signal, higher values means less noise
y_obs <- t(apply(y_prob, 1, function(p) rmultinom(1, clean, p) / clean))

# make sure we have a max in each row for each group
table(apply(y_obs, 1, function(x) which(x == max(x))[1]))

# make sure all the values are close to 1
print(all.equal(rowSums(y_prob), rep(1, N)))
print(all.equal(rowSums(y_obs), rep(1, N)))

# function for predicting
ratio2Pred <- function(X, model, draws=NULL){
    G <- length(model$par) / ncol(X) + 1
    if(!is.null(draws)){
        # recursive process for simulations
        betas_ <- rmvnorm(draws, model$par, solve(model$hessian))
        y_prob <- array(0, dim=c(nrow(X), G, draws))
        for(d in 1:draws){
            m <- list(par=betas_[d,]) # treat the sims like true values
            y_prob[,,d] <- ratio2Pred(X, m) # rerun the function with sim
        }
    }
    else{
        # transform vector of beta paramters into matrix for ease of use
        betas_ <- matrix(model$par, ncol=G-1, nrow=ncol(X))
        # get the ratios for the groups vs ref
        ratio_ref <- exp(X %*% betas_)
        # the response var is transformed to probability space
        p1 <- apply(ratio_ref, 1, function(x) 1 / (1 + sum(x)))
        y_prob <- cbind(p1, ratio_ref * p1)
    }
    return(y_prob)
}

# get the estimates
optimMultiVec <- function(params, X, Y){
    # apply the log likelihood to the observations
    -sum(log(ratio2Pred(X, list(par=params))) * Y)
}

fitMultiVec <- function(X, Y){
    # problematic if not all groups are well defined
    # probably can get around this with sampling
    y <- apply(y_obs, 1, function(x) which(x == max(x))[1])

    # get starting values from regular ass multinom
    initDF <- as.data.frame(cbind(y=y, X))
    pstart <- c(coef(multinom(y ~ . + 0, data=initDF, trace=FALSE)))

    # optimizeeeee
    optim(pstart, optimMultiVec, X=X, Y=Y, hessian=TRUE, method="BFGS")
}

# test em
testFit <- fitMultiVec(covs, y_obs)

testFit$convergence # make sure we get 0 exit status
# visually inspect betas not bad!
betas
matrix(testFit$par, nrow = M+1)

# we can also get draws of probability predictions
dim(ratio2Pred(covs, testFit, draws = 1000))
