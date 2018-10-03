rm(list=ls())
pacman::p_load(INLA, ggplot2, data.table, lattice, arm, dplyr, TMB, ar.matrix)
set.seed(124)
# compare the INLA Q matrix vs the by hand to make sure we are on the same page

mesh2DF <- function(x){
    M <- length(proj$x)
    DT <- data.table(x=rep(proj$x, M), y=rep(proj$y, each=M), 
                     obs=c(inla.mesh.project(proj, field=x)))
    DT
}

n <- 1000
loc <- matrix(runif(n*2), n, 2)
mesh <- inla.mesh.create(loc, refine=list(max.edge=0.05))

plot(mesh)
points(loc[,1], loc[,2], col="red", pch=20)
proj <- inla.mesh.projector(mesh, dims=c(500, 500))

sigma0 <-  .2   ## Standard deviation
range0 <- .1 ## Spatial range
kappa0 <- sqrt(8) / range0
tau0 <- 1/(sqrt(4*pi)*kappa0*sigma0)
spde <- inla.spde2.matern(mesh)

##############
Q2 <- tau0**2 * (kappa0**4 * spde$param.inla$M0 + 
                     2 * kappa0**2 *spde$param.inla$M1 + spde$param.inla$M2)
x_ <- as.vector(sim.AR(1, Q2))
x <- x_ - mean(x_)

mesh2DF(x) %>%
    ggplot(aes(x, y, z=obs)) +
    geom_raster(aes(fill = obs)) + 
    theme_void() + 
    lims(y=c(0,1), x=c(0,1)) + 
    scale_fill_distiller(palette = "Spectral")

cov1 <- rnorm(n)

p <- invlogit(-1 + .2 * cov1 + x[mesh$idx$loc])
denom <- rpois(n, 100)

hist(p)
hist(denom)

DT <- data.table(
    y=rbinom(n, denom, p),
    cov=cov1, id=1:n, denom=denom,
    geo=mesh$idx$loc-1,
    lon=loc[,1], lat=loc[,2])
summary(DT)

spatform <- y ~ cov1 + f(id, model=spde)

system.time(spatmodel <- inla(
    spatform, 
    data=DT,
    control.compute=list(config=TRUE),
    family="binomial",
    Ntrials = DT$denom))

summary(spatmodel)
spatialhat <- inla.spde2.result(spatmodel, "id", spde)
xhat <- spatialhat$summary.values$`0.5quant`
bind_rows(
    mutate(mesh2DF(x), state="Observed"),
    mutate(mesh2DF(xhat), state="Estimated INLA")) %>%
    ggplot(aes(x, y, z=obs)) +
    geom_raster(aes(fill = obs)) +
    theme_void() +
    lims(y=c(0,1), x=c(0,1)) +
    scale_fill_distiller(palette = "Spectral") +
    facet_wrap(~state)

mesh2DF(spatialhat$summary.values$sd) %>%
    ggplot(aes(x, y, z=obs)) +
    geom_raster(aes(fill = obs)) +
    theme_void() +
    scale_fill_distiller(palette = "Spectral") +
    geom_point(aes(x=lon, y=lat, z=NULL, fill=NULL), data=DT)

exp(spatmodel$summary.hyperpar$mean[1:2])
c(tau0, kappa0)

# now run the TMB model using ./st.cpp
setwd("~/Documents/re_simulations/ppp/")

runModel <- function(DT, recompile=FALSE, symbolic=TRUE, draws=1000){
    Data <- list(
        y=DT$y, x=DT$cov, geo=DT$geo, denom=DT$denom,
        M0=spde$param.inla$M0, M1=spde$param.inla$M1, M2=spde$param.inla$M2)
    
    Params <- list(
        beta0=0, beta1=0, log_tau=0, log_kappa=0, z=rep(0, nrow(Q2))
    )
    
    # compile the code if not there
    model <- "ppp"
    if(recompile){
        if (file.exists(paste0(model, ".so"))) file.remove(paste0(model, ".so"))
        if (file.exists(paste0(model, ".o"))) file.remove(paste0(model, ".o"))
        if (file.exists(paste0(model, ".dll"))) file.remove(paste0(model, ".dll"))
    }
    compile(paste0(model, ".cpp"))
    
    dyn.load(dynlib(model))
    Obj <- MakeADFun(data=Data, parameters=Params, DLL=model, random="z")
    if(symbolic){
        runSymbolicAnalysis(Obj)
    }
    Opt <- nlminb(
        start=Obj$par,
        objective=Obj$fn,
        gradient=Obj$gr,
        control=list(eval.max=1e4, iter.max=1e4))
    sdrep <- sdreport(Obj, getJointPrecision=TRUE)
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
    list(obj=Obj, opt=Opt, sd=sdrep, z=zDF)
}

system.time(tmbModel <- runModel(DT))

bind_rows(
    mutate(mesh2DF(x), state="Observed"),
    mutate(mesh2DF(tmbModel$z$mu), state="Estimated TMB"),
    mutate(mesh2DF(xhat), state="Estimated INLA")) %>%
    ggplot(aes(x, y, z=obs)) +
    geom_raster(aes(fill = obs)) +
    theme_void() +
    lims(y=c(0,1), x=c(0,1)) +
    scale_fill_distiller(palette = "Spectral") +
    facet_wrap(~state, nrow=2)

mesh2DF(tmbModel$z$sd) %>%
    ggplot(aes(x, y, z=obs)) +
    geom_raster(aes(fill = obs)) +
    theme_void() +
    scale_fill_distiller(palette = "Spectral") +
    geom_point(aes(x=lon, y=lat, z=NULL, fill=NULL), data=DT)
