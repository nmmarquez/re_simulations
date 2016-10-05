set.seed(123)
# library to run multinomial regression
library(nnet)
rm(list=ls())

N <- 100000 # number of observations
M <- 2 # number of non-intercept covariates
G <- 5 # number of groups

# generate the covariates
covs <- sapply(1:M, function(x) rnorm(N))
covs <- cbind(rep(1, N), covs)

# generate the betas
betas <- sapply(1:(G-1), function(x) runif(M+1, min = -.4, max=.4))

# the response var in log_ratio_space
ratio_ref <- exp(covs %*% betas)

# the response var in probability space
p1 <- apply(ratio_ref, 1, function(x) 1 / (1 + sum(x)))
y_prob <- cbind(p1, ratio_ref * p1)

# make sure all the values are close to 1
print(all.equal(rowSums(y_prob), rep(1, N)))

# sample once from row with weighted probabilities
y <- apply(y_prob, 1, function(x) sample(1:G, 1, prob=x))
table(y)

# put info in a data frame for model function conveinience 
df <- as.data.frame(cbind(y, covs[,2:(M+1)]))
df$y2 <- relevel(as.factor(df$y), "1")

# build expression and run model
expr <- paste0("y2 ~ ", paste("V", 2:(M+1), collapse= " + ", sep=""))
print("Run time of nnet model for multinomial log-linear regression...")
system.time(mnlr <- multinom(formula(expr), data=df, trace=F))

# look at the betas
print("rmse of betas when using nnet")
mean((t(betas) - summary(mnlr)$coefficients)**2)

# run the model using TMB
library(TMB)
setwd("~/Documents/tmp/mnlr/")

# compile the cpp
model_name <- "mnlr"
#if (file.exists(paste0(model_name, ".so"))) file.remove(paste0(model_name, ".so"))
#if (file.exists(paste0(model_name, ".o"))) file.remove(paste0(model_name, ".o"))
#if (file.exists(paste0(model_name, ".dll"))) file.remove(paste0(model_name, ".dll"))
compile(paste0(model_name, ".cpp"))

# load the model
dyn.load(dynlib(model_name))

# set the starting parameter points and data
Params <- list(betas=matrix(0, nrow=M+1, ncol=G-1))
Data <- list(group=y-1, covs=covs)

# build and optimize the objective function
Obj <- MakeADFun(data=Data, parameters=Params, DLL=model_name, silent=T)
Obj$env$tracemgc <- FALSE
Obj$env$inner.control$trace <- FALSE
Obj$env$silent <- TRUE
print("Run time of TMB model for multinomial log-linear regression...")
system.time(Opt <- nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr))

# check out the betas actually a bit better than the nnet implementation
# although a lot slower
Report <- Obj$report()
print("rmse of betas when using TMB")
mean((betas - Report$betas)**2)