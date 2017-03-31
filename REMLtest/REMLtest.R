# load your packages
rm(list=ls())
pacman::p_load(data.table, lme4, TMB)

setwd("~/Documents/re_simulations/REMLtest/")

# create a groups variable that TMB can use 
iris$Group <- as.numeric(as.factor(iris$Species)) - 1

# REML and ML in lme4 and check the diff in output
REMLlm <- lmer(Sepal.Length ~ Sepal.Width + (1|Group), data=iris)
summary(REMLlm)
MLlm <- lmer(Sepal.Length ~ Sepal.Width + (1|Group), data=iris, REML=FALSE)
summary(MLlm)
# both of these are really different than not inlcuding the random effect
# this is pretty normal if you consider something like simpsons paradox
reglm <- lm(Sepal.Length ~ Sepal.Width, data=iris)
summary(reglm)

# lets emulate this in tmb
model <- "REMLtest"
if (file.exists(paste0(model, ".so"))) file.remove(paste0(model, ".so"))
if (file.exists(paste0(model, ".o"))) file.remove(paste0(model, ".o"))
if (file.exists(paste0(model, ".dll"))) file.remove(paste0(model, ".dll"))
compile(paste0(model, ".cpp"))

run_model <- function(REML=TRUE, model_name=model){
    dyn.load(dynlib(model_name))
    
    # set the starting parameter points and data
    Params <- list(b0=0, b1=0, log_sigma=0, log_re_sigma=0, 
                   zeta=rep(0, length(unique(iris$Group))))
    Data <- list(y=iris$Sepal.Length, x=iris$Sepal.Width, group=iris$Group)
    
    # build and optimize the objective function
    ran <- c("zeta")
    if(REML){
        # if REML the betas are treated as "random effects" really their value 
        # is just treated as unknown and coming from a random variable so 
        # we integrate it out of the likelihood with our other random variables
        ran <- c(ran, "b0", "b1")
    }
    Obj <- MakeADFun(data=Data, parameters=Params, DLL=model_name, silent=T, 
                     random=ran)
    Obj$env$tracemgc <- FALSE
    Obj$env$inner.control$trace <- FALSE
    Obj$env$silent <- TRUE
    print("Run time of TMB model for multinomial log-linear regression...")
    system.time(Opt <- nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr))
    Opt$convergence
    
    Report <- Obj$report()
    dyn.unload(dynlib(model_name))
    Report
}

# set a pretty good tolerance
tol <- 10**-6

# TMB REML matches with LMER REML
REMLtmb <- run_model()
all.equal(REMLtmb$b0, REMLlm@beta[1], tolerance=tol)
all.equal(REMLtmb$b1, REMLlm@beta[2], tolerance=tol)

# same with non reml version
MLtmb <- run_model(REML=FALSE)
all.equal(MLtmb$b0, MLlm@beta[1], tolerance=tol)
all.equal(MLtmb$b1, MLlm@beta[2], tolerance=tol)
