rm(list=ls())
library(TMB)
library(tidyverse)
library(ar.matrix) # devtools::install_github("nmmarquez/ar.matrix")

# simulate some easy peasy data
N <- 500
sigma_epsilon <- .2
sigma_ar <- .5
rho_ar <- .92
beta0 <- 5

# simulate with seed
set.seed(123)
yobs <- rnorm(N, beta0 + c(r.AR1(1, N, sigma_ar, rho_ar)), sigma_epsilon)

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

models <- c(sequence=0, GMRF=1, builtAR1=2, MVNM=3)

modelFits <- lapply(models[1:4], function(i) runModel(yobs, i))

sapply(modelFits, function(x) x$opt$par)
sapply(modelFits, function(x) x$runtime)
