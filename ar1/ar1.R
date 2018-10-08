rm(list=ls())
library(TMB)
library(tidyverse)


# simulate some easy peasy data
N <- 50000
sigma_epsilon <- .2
sigma_ar <- .5
rho_ar <- .92
beta0 <- 5

# simulate with seed
set.seed(123)
handAR <- function(N, sigma_ar, rho_ar){
    vec <- rep(0, N)
    vec[1] <- rnorm(1, 0, sqrt(sigma_ar^2 / (1 - rho_ar^2)))
    for(j in 2:N){
        vec[j] <- rnorm(1, rho_ar*vec[j-1], sigma_ar)
    }
    vec
}

yobs <- rnorm(N, beta0 + handAR(N, sigma_ar, rho_ar), sigma_epsilon)

# compile if needed
compile("ar1.cpp")

## Model Run Function

runModel <- function(y, option){
    if(!(option %in% 0:3)){
        stop("Option must be 0, 1, 2, or 3")
    }
    
    dyn.load(dynlib("ar1"))
    
    start_time <- Sys.time()
    
    Obj <- MakeADFun(
        data = list(yobs = y, option = option),
        parameters = list(
            log_sigma_epsilon = 0,
            log_sigma_ar = 0,
            logit_rho_ar = 0,
            beta = 0,
            zeta = rep(0, N)
        ),
        random = "zeta",
        DLL = "ar1"
    )
    
    Opt <- nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr)
    
    sdrep <- sdreport(Obj, getJointPrecision = T)
    
    run_time <- Sys.time() - start_time
    
    dyn.unload(dynlib("ar1"))
    
    list(obj=Obj, opt=Opt, sdrep=sdrep, runtime=run_time)
}

models <- c(sequence=0, GMRF=1, builtAR1=2)#, MVNM=3)

modelFits <- lapply(models[1:4], function(i) runModel(yobs, i))

sapply(modelFits, function(x) x$opt$par)
sapply(modelFits, function(x) x$runtime)
