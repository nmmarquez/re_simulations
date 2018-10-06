# We want to get an idea if its feasible to constrain values in TMB, especially
# random effects, with a sum to zero constraint. Here we test the sum to zero 
# constraint for a set of unstructured random effects with a simple linear
# regression. The regression containes two fixed effects and a number of random 
# effects equal to 10% of the sample size.

rm(list=ls())
library(TMB)
library(lme4)
library(parallel)
library(dplyr)

args <- commandArgs(trailingOnly=TRUE)
print(args)
N <- as.integer(args[1])
M <- N/10
const <- as.integer(args[2])
seed <- as.integer(args[3])

setwd("~/Documents/re_simulations/constrainTMB/")

runSimulation <- function(N, M, const, a=3, b=5, sigmaObs=2, sigmaZ=.7, s=123){
    set.seed(s)
    x <- runif(N) # Simulate the covariate values
    
    groups <- sample(1:M, N, replace = T) # randomly sample groups for each obs
    z <- rnorm(M, sd=sigmaZ) # simulate the random effects
    y <- rnorm(N, a + x * b + z[groups], sigmaObs) # simulate the obs
    
    dyn.load(dynlib("constrain")) # load tmb library
    
    time_1 <- system.time(Obj <- MakeADFun(
        data=list(Y=y, x=x, group=groups-1, constrain=const),
        parameters=list(a=0, b=0, logSigma=0, logSigmaZ=0, 
                        z=ifelse(const==1, M-1, M)),
        random="z",
        DLL="constrain"))
    
    time_2 <- system.time(Opt <- nlminb(
        start=Obj$par, objective=Obj$fn, 
        gradient=Obj$gr,
        control=list(eval.max=1e4, iter.max=1e4)))
    
    time_3 <- system.time(sdrep <- sdreport(Obj, getJointPrecision = T))
    
    dyn.unload(dynlib("constrain"))
    
    time_all <- time_1 + time_2 + time_3
    
    tibble(
        b0=Opt$par[1], b1=Opt$par[2],
        sigma.resid=exp(Opt$par[3]), sigma.random=exp(Opt$par[4]),
        time.user=time_all[1], time.system=time_all[2],
        time.elapsed=time_all[3], n=N, m=M, constrain=const)
}

simRezDF <-runSimulation(N, M, const, s=seed)

write.csv(
    simRezDF,
    paste0(
        "~/Documents/re_simulations/constrainTMB/results/", 
        N, "_", const, "_", seed, ".csv"),
    row.names = F
)
