rm(list=ls())

library(boot)
library(tidyverse)
library(sf)
library(sp)
library(INLA)
library(ar.matrix)
library(TMB)

mesh2DF <- function(model){
    if(class(model) == "numeric"){
        x <- model
        sdx <- rep(NA, length(x))
    }
    else{
        x <- model$z$mu
        sdx <- model$z$sd
    }
    M <- length(proj$x)
    DT <- data.frame(x=rep(proj$x, M), y=rep(proj$y, each=M), 
                     obs=c(inla.mesh.project(proj, field=x)),
                     sdx=c(inla.mesh.project(proj, field=sdx)))
    expand.grid(proj$x, proj$y) %>%
        SpatialPoints(US.df@proj4string) %>%
        over(US.df) %>%
        cbind(DT) %>%
        mutate(aproj=1:n(), obsField=!is.na(id)) %>%
        group_by(obsField) %>%
        mutate(obsAproj=1:n()) %>%
        ungroup %>%
        as_tibble %>%
        mutate(stateID=group_indices(., STATEFP)) %>%
        mutate(stateID=ifelse(!obsField, NA, stateID)) %>%
        mutate(obsAproj=ifelse(!obsField, NA, obsAproj))
}

# what is the shape that we are dealing with?
row.names(US.df@data) <- 1:nrow(US.df)
US.df$id <- 1:nrow(US.df)
USDF <- fortify(US.df, region="id") %>%
    mutate(id=as.numeric(id)) %>%
    left_join(US.df@data, by="id")

set.seed(123)
randomSPDF <- spsample(US.df, 800, "random")
randomSPDF$long <- randomSPDF@coords[,1]
randomSPDF$lat <- randomSPDF@coords[,2]

plot(mesh <- inla.mesh.2d(
    randomSPDF, 
    cutoff=.5,
    max.edge=c(50, 500)))
proj <- inla.mesh.projector(mesh, dims=c(250, 250))

beta0 <- -1
sigma0 <-  .6   ## Standard deviation
range0 <- 1.5 ## Spatial range
kappa0 <- sqrt(8) / range0
tau0 <- 1/(sqrt(4*pi)*kappa0*sigma0)
spde <- inla.spde2.matern(mesh)
m <- 9 # time periods
rho <- .93 # temporal autocorrelation
Qspde <- tau0**2 * (kappa0**4 * spde$param.inla$M0 + 
                        2 * kappa0**2 *spde$param.inla$M1 + spde$param.inla$M2)

Q <- kronecker(Q.AR1(m, 1, rho), Qspde)
x_ <- matrix(data=c(sim.AR(1, Q)), nrow=mesh$n, ncol=m)
x <- x_ - mean(x_)

fieldPlot <- bind_rows(lapply(1:m, function(j){
    mesh2DF(x[,j]) %>%
        mutate(time=j)})) %>%
    filter(obsField) %>%
    mutate(p=inv.logit(beta0 + obs)) %>%
    ggplot(aes(x, y, z=p)) +
    geom_raster(aes(fill=p)) + 
    coord_equal() +
    theme_void() + 
    scale_fill_distiller(palette = "Spectral") +
    facet_wrap(~time)

fieldPlot

## same plot with boundaries added
fieldPlot +
    geom_path(
        aes(long,lat, group=group, fill=NULL, z=NULL),
        color="black",
        size=.1,
        data=USDF)

N <- 2000
obsPoints <- spsample(US.df, N, "random")
obsPoints$time <- sample(c(1:4, 7:9), N, replace=T)
AprojPoint <- inla.spde.make.A(
    mesh = mesh,
    loc = obsPoints,
    group = obsPoints$time)

obsDF <- tibble(long=obsPoints@coords[,1], lat=obsPoints@coords[,2]) %>%
    mutate(re=as.vector(AprojPoint %*% c(x)), denom=rpois(N, 100)) %>%
    mutate(prob=inv.logit(beta0 + re), obs=rbinom(N, denom, prob))

runModel <- function(DFpoint=NULL, recompile=F, verbose=F, draws=1000){
    model <- "st"
    compile(paste0(model, ".cpp"))
    Data <- list(
        yPoint=DFpoint$obs, denomPoint=DFpoint$denom, AprojPoint=AprojPoint,
        M0=spde$param.inla$M0, M1=spde$param.inla$M1, M2=spde$param.inla$M2,
        timem=m)
    
    Params <- list(
        beta0=0, log_tau=0, log_kappa=0, z=rep(0, nrow(Q)), logit_rho=0
    )
    
    dyn.load(dynlib(model))
    startTime <- Sys.time()
    Obj <- MakeADFun(data=Data, parameters=Params, DLL=model, random="z",
                     silent=!verbose)
    Obj$env$tracemgc <- verbose
    Obj$env$inner.control$trace <- verbose
    symbolic <- T
    if(symbolic){
        nah <- capture.output(runSymbolicAnalysis(Obj))
    }
    Opt <- nlminb(
        start=Obj$par,
        objective=Obj$fn,
        gradient=Obj$gr,
        control=list(eval.max=1e4, iter.max=1e4))
    sdrep <- sdreport(Obj, getJointPrecision=TRUE)
    runtime <- Sys.time() - startTime
    if(attr(runtime, "units") == "mins"){
        runtime <- as.numeric(runtime) * 60
    }
    else if(attr(runtime, "units") == "secs"){
        runtime <- as.numeric(runtime)
    }
    zindex <- "z" == row.names(sdrep$jointPrecision)
    Qz <- sdrep$jointPrecision[zindex,zindex]
    Zdraws <- sim.AR(draws, Qz)
    Report <- Obj$report()
    zDF <- tibble(
        mu=Report$z,
        sd=apply(Zdraws, 2, sd),
        lwr=apply(Zdraws, 2, quantile, probs=.025),
        upr=apply(Zdraws, 2, quantile, probs=.975)
    )
    return(list(obj=Obj, opt=Opt, z=zDF, runtime=runtime, sd=sdrep))
}

modelRES <- runModel(obsDF, recompile=FALSE, verbose=TRUE)
