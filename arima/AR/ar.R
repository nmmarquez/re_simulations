rm(list = ls())
pacman::p_load(TMB)

# test out the moving average stuff
N <- 10000
mu <- 2
p <- c(1.)
Np <- length(p)
x <- rep(0, N)

for(i in 1:Np){
    if(p != 1){
        x[i] <- mu / sum(c(1, -1 * p))
    }
    else{
        x[i] <- mu * i
    }
}

for(i in (Np + 1):N){
    x[i] <- mu + sum(x[(i-1):(i-Np)] * p) + rnorm(1, sd=1)
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
Data <- list(x=x)

# build and optimize the objective function
Obj <- MakeADFun(data=Data, parameters=Params, DLL=model_name)
system.time(Opt <- nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr))
Opt$convergence

(sdrep <- sdreport(Obj))
rep <- Obj$report()
arima(x, c(1,0,0))

# differences are small
arima(x, c(Np,0,0))$coef[paste0("ar", 1:Np)] - sdrep$value[-1]

