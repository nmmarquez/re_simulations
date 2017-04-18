rm(list = ls())
pacman::p_load(TMB)

# test out the moving average stuff
N <- 10000
w <- rnorm(N)
mu <- 0
q <- c(.8)
Nq <- length(q)
x <- rep(0, N)

for(i in (Nq + 1):N){
    x[i] <- mu + w[i] + sum(w[(i-Nq):(i-1)] * q)
}
plot(x, type="l")

mean(x)
var(x)
acf(x)

setwd("~/Documents/re_simulations/arima/MA/")

model_name <- "ma"
if (file.exists(paste0(model_name, ".so"))) file.remove(paste0(model_name, ".so"))
if (file.exists(paste0(model_name, ".o"))) file.remove(paste0(model_name, ".o"))
if (file.exists(paste0(model_name, ".dll"))) file.remove(paste0(model_name, ".dll"))
compile(paste0(model_name, ".cpp"))

# load the model
dyn.load(dynlib(model_name))

# set the starting parameter points and data
Params <- list(w=rep(0, N), log_sigma=0, log_sigmaw=0, q=0, mu=0)
Data <- list(x=x)

arima(x, c(0,0,1))

# build and optimize the objective function
Obj <- MakeADFun(data=Data, parameters=Params, DLL=model_name, random="w")
system.time(Opt <- nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr))
Opt$convergence

(sdrep <- sdreport(Obj))
rep <- Obj$report()

