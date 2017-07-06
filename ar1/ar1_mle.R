rm(list=ls())
set.seed(123)
# load the packages to evaluate multivariate normals with covariance and 
# precision matrix
pacman::p_load(mvtnorm, LaplacesDemon)

# simulate data using R's built in ARIMA simulation
X <- as.vector(arima.sim(n=100, list(ar=0.8897), sd=0.2))

# function to build a precision matrix of AR1 process given rho and sigma
QAR <- function(N, rho, sigma){
    Q <- matrix(0, nrow=N, ncol=N)
    Q[1, 1] <- 1. / sigma**2
    for(i in 2:N){
        Q[i, i] <- (1. + rho**2) / sigma**2
        Q[i-1, i] <- (-1. * rho) / sigma**2
        Q[i, i-1] <- (-1. * rho) / sigma**2
    }
    Q[N, N] <- 1. / sigma**2
    return(Q)
}

# function to build a covariance matrix of AR1 process given rho and sigma
SigmaAR <- function(N, rho, sigma){
    H <- abs(outer(1:N, 1:N, "-"))
    return(sigma**2 * rho**H / (1 - rho**2))
}

densf <- function(pars, x=X, option=1){
    rho <-  1. / (1. + exp(-1. * pars[1]))
    sigma <- exp(pars[2])
    N <- length(x)
    nll <- 0.
    # option 1 uses a for loop to estimate the AR1 process
    if(option == 1){
        for(i in 2:N){
            nll <- nll - dnorm(x[i], rho * x[i-1], sigma, log=TRUE)
        }
    }
    # option 2 uses a covariance matrix and the dmvnorm function from mvtnorm
    else if(option == 2){
        Sigma <- SigmaAR(N, rho, sigma)
        nll <- -1 * dmvnorm(x, sigma=Sigma, log=TRUE)
    }
    # option 3 uses a precision matrix and the dmvn function from LaplacesDemon
    else if(option == 3){
        Q <- QAR(N, rho, sigma)
        nll <- -1 * dmvnp(x, rep(0, N), Q, log=TRUE)
    }
    return(nll)
}

# utility function to transform parameters
transform_pars <- function(pars){
    c(rho=1. / (1. + exp(-1. * pars[1])), sigma=exp(pars[2]))
}

transform_pars(optim(c(0,0), densf, option=1)$par) # loop methods
transform_pars(optim(c(0,0), densf, option=2)$par) # covariance method
transform_pars(optim(c(0,0), densf, option=3)$par) # precision method
can.arima <- arima(X, c(1,0,0), include.mean=FALSE) # base r arima method
c(rho=can.arima$coef[[1]], sigma=can.arima$sigma2**.5)
