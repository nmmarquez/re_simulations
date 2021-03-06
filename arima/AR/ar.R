rm(list = ls())
pacman::p_load(TMB)

# test out the moving average stuff
N <- 1000
mu <- 0
p <- c(.8, -.3)
Np <- length(p)
x <- rep(0, N)
d <- 2

for(i in 1:Np){
    if(abs(sum(p)) < 1){
        x[i] <- mu / sum(c(1, -1 * p))
    }
    else{
        x[i] <- mu * i
    }
}


for(i in (Np + 1):N){
    x[i] <- mu + sum(x[(i-1):(i-Np)] * p) + rnorm(1, sd=1)
}

if(d > 0){
    for(i in 1:d){
        x <- c(0, cumsum(x))
    }
}

plot(x, type="l")

setwd("~/Documents/re_simulations/arima/AR/")

model_name <- "ar"
# if (file.exists(paste0(model_name, ".so"))) file.remove(paste0(model_name, ".so"))
# if (file.exists(paste0(model_name, ".o"))) file.remove(paste0(model_name, ".o"))
# if (file.exists(paste0(model_name, ".dll"))) file.remove(paste0(model_name, ".dll"))
compile(paste0(model_name, ".cpp"))

# load the model
dyn.load(dynlib(model_name))

# set the starting parameter points and data
Params <- list(mu=0, log_sigma=0, p=rep(0, Np))
Map <- list()
x_ <- x
if(d > 0){
    x_ <- diff(x, differences=d)
    Map[["mu"]] <- as.factor(NA)
}
Data <- list(x=x_)

# build and optimize the objective function
Obj <- MakeADFun(data=Data, parameters=Params, DLL=model_name, map=Map)
system.time(Opt <- nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr))
Opt$convergence

(sdrep <- sdreport(Obj, getJointPrecision = T))
rep <- Obj$report()

# differences are small
arima(x, c(Np,d,0))$coef[paste0("ar", 1:Np)] - sdrep$value[-1]
arima(x, c(Np,d,0))$coef

# estimated mean value from model
sdrep$value["mu"] / (1 - sum(sdrep$value[-1]))


# # RW with drift?
# source("~/Documents/re_simulations/arima/AR/snipy/utilities.R")
# set.seed(1234)
# DF <- sim_ar_data(N=1000, p= 1., mu = 2)
# plot(DF$obs)
# forecast::auto.arima(DF$obs)
# forecast::Arima(DF$obs, c(0,1,0), include.drift=T)
# arima(DF$obs, c(1,0,0), xreg=1:nrow(DF), include.mean=F)
# run_arima_TMB(DF$obs, 1)
# 
# testlm(DF$obs ~ DF$time)
