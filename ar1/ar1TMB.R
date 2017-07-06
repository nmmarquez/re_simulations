rm(list=ls())
pacman::p_load(TMB, data.table)

X <- as.vector(arima.sim(n=100, list(ar=0.8897), sd=0.2))

# now run the TMB model using ./st.cpp
setwd("~/Documents/re_simulations/ar1/")

model <- "ar1TMB"
if (file.exists(paste0(model, ".so"))) file.remove(paste0(model, ".so"))
if (file.exists(paste0(model, ".o"))) file.remove(paste0(model, ".o"))
if (file.exists(paste0(model, ".dll"))) file.remove(paste0(model, ".dll"))
compile(paste0(model, ".cpp"))

run_model <- function(option, X_=X){
    # set the data
    Data <- list(yobs=X_, option=option)
    
    # set the param starting points
    Params <- list(logsigma=0, logitrho=0)
    
    # load and optimize
    dyn.load(dynlib(model))
    Obj <- MakeADFun(data=Data, parameters=Params, DLL=model, silent=T)
    Obj$env$tracemgc <- F
    Obj$env$inner.control$trace <- F
    Opt <- nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr)
    Report <- Obj$report()
    dyn.unload(dynlib(model))
    Report
}


ar1model <- arima(X, order=c(1,0,0), include.mean=F)

rbindlist(list(
    data.frame(sigma=ar1model$sigma2**.5, rho=ar1model$coef, method=0),
    data.frame(c(run_model(1)), list(method=1)),
    data.frame(c(run_model(2)), list(method=2)),
    data.frame(c(run_model(3)), list(method=3)),
    data.frame(c(run_model(4)), list(method=4))
))

