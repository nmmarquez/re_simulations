rm(list=ls())
set.seed(123)
# run the model using TMB
library(TMB)
library(MASS)
setwd("~/Documents/re_simulations/nb_compare/")

reload_model <- function(model_name){
    file_exts <- c(".so", ".o", ".dll")
    for (x in file_exts){
        if (file.exists(paste0(model_name, x))){
            file.remove(paste0(model_name, x))
        }
    }
    compile(paste0(model_name, ".cpp"))
}

# NB1 model which is a mixture model of poisson with extra gamma distributed
# paramter alpha
model_name <- "nb1"
reload_model(model_name)

N <- 10000 # sample size
alpha <- .03 # overdispersion parameter
betas <- c(2, .5) # betas
M <- length(betas) # number of covariates
nu <- log(rgamma(N, shape=1/alpha, scale=alpha)) # overdispersion realizations
X <- cbind(rep(1, N), sapply(2:M, function(x) rnorm(N))) # sim covariates
mu <- exp((X %*% betas) + nu) # the linear model with a log link function
y <- rpois(N, mu) # simulate the data from a poisson distribution
hist(y)

# we can run this as a regular poisson model if we map nu to NA in the objective
Map <- list()
Map[["nu"]] <- factor(array(NA, dim=c(N)))

# load the model
dyn.load(dynlib(model_name))

# set the starting parameter points and data
Params <- list(betas=rep(0, M), log_alpha=log(.7), nu=rep(0, N))
Data <- list(X=X, y=y)

# build and optimize the objective function
Obj <- MakeADFun(data=Data, parameters=Params, random="nu", DLL=model_name)
# We make nu random here as our unobserved random variable 
Opt <- nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr)
# make sure we got that zero exit status
Opt$convergence

# check our results
Report <- Obj$report()
Report$betas
Report$alpha

# get the joint precision matrix of the latent random effects
(sdrep <- sdreport(Obj, getJointPrecision = T))

dyn.unload(dynlib(model_name))

# Now we will simulate for teh NB2 model using R's built in rndbinom function
model_name <- "nb2"
reload_model(model_name)

theta <- 1 / alpha # overdispersion parameter which is the inverse of alpha
y2 <- rnbinom(N, mu=exp(X %*% betas), size=theta) # we can use teh same covs
hist(y2)

# load the model
dyn.load(dynlib(model_name))

# set the starting parameter points and data
Params <- list(betas=rep(0, M), log_theta=log(1/1))
Data <- list(X=X, y=y2)

# build and optimize the objective function
Obj <- MakeADFun(data=Data, parameters=Params, DLL=model_name)
Opt <- nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr)
Opt$convergence

Report <- Obj$report()
Report$betas
Report$theta
