rm(list=ls())
pacman::p_load(INLA, ggplot2, data.table, lattice, arm, dplyr, TMB, 
               ar.matrix, dtplyr)
set.seed(124)

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
    DT <- data.table(x=rep(proj$x, M), y=rep(proj$y, each=M), 
                     obs=c(inla.mesh.project(proj, field=x)),
                     sdx=c(inla.mesh.project(proj, field=sdx)))
    expand.grid(proj$x, proj$y) %>%
        SpatialPoints(US.df@proj4string) %>%
        over(US.df) %>%
        as.data.table %>%
        cbind(DT) %>%
        mutate(aproj=1:n(), obsField=!is.na(id)) %>%
        group_by(obsField) %>%
        mutate(obsAproj=1:n()) %>%
        ungroup %>%
        as_tibble %>%
        mutate(stateID=group_indices(., STATEFP)) %>%
        mutate(stateID=ifelse(!obsField, NA, stateID)) %>%
        mutate(obsAproj=ifelse(!obsField, NA, obsAproj)) %>%
        as.data.table
}

# what is the shape that we are dealing with?
row.names(US.df@data) <- 1:nrow(US.df)
US.df$id <- 1:nrow(US.df)
USDF <- fortify(US.df, region="id") %>%
    mutate(id=as.numeric(id)) %>%
    left_join(US.df@data, by="id")

randomSPDF <- spsample(US.df, 800, "random")
randomSPDF$long <- randomSPDF@coords[,1]
randomSPDF$lat <- randomSPDF@coords[,2]

plot(mesh <- inla.mesh.2d(
    randomSPDF, 
    cutoff=.1,
    max.edge=c(50, 500)))
proj <- inla.mesh.projector(mesh, dims=c(500, 500))

beta0 <- -1
sigma0 <-  .6   ## Standard deviation
range0 <- 1.5 ## Spatial range
kappa0 <- sqrt(8) / range0
tau0 <- 1/(sqrt(4*pi)*kappa0*sigma0)
spde <- inla.spde2.matern(mesh)

Q <- tau0**2 * (kappa0**4 * spde$param.inla$M0 + 
                     2 * kappa0**2 *spde$param.inla$M1 + spde$param.inla$M2)
x_ <- as.vector(sim.AR(n=1, Q))
x <- x_ - mean(x_)

simValues <- mesh2DF(x)

simValues %>%
    filter(obsField) %>%
    mutate(p=invlogit(beta0 + obs)) %>%
    ggplot(aes(x, y, z=p)) +
    geom_raster(aes(fill=p)) + 
    coord_equal() +
    theme_void() + 
    scale_fill_distiller(palette = "Spectral") +
    geom_path(
        aes(long,lat, group=group, fill=NULL, z=NULL),
        color="black",
        size=.1,
        data=USDF)

# assume equal pop weight across pixels
countyRiskDF <- simValues %>%
    filter(obsField) %>%
    group_by(id) %>%
    summarize(p=mean(invlogit(beta0 + obs), na.rm=T), cells=n())

countyRiskDF %>%
    right_join(USDF, by="id", copy=T) %>%
    ggplot +
    aes(long,lat,group=group,fill=p) + 
    geom_polygon() +
    geom_path(color="black", size=.1) +
    coord_equal() + 
    theme_void() +
    scale_fill_distiller(palette = "Spectral")

# how does the estimation at multiple levels very our probability estimate
# distribution?
idtest <- 392
simArea <- function(N, M, id=idtest){
    obsPoints <- spsample(US.df[US.df$id==id,], N, "random")@coords
    AprojPoint <- inla.spde.make.A(mesh=mesh, loc=obsPoints)
    obsDF <- data.table(long=obsPoints[,1], lat=obsPoints[,2]) %>%
        mutate(re=as.vector(AprojPoint %*% x), denom=M) %>%
        mutate(prob=invlogit(beta0 + re), obs=rbinom(N, denom, prob))
    pCounty <- countyRiskDF$p[countyRiskDF$id==id]
    c(mean(obsDF$obs / M), rbinom(1, N*M, pCounty)/(N*M))
}

# nSims <- 1000
# data.frame(
#     prob=c(sapply(1:nSims, function(i) simArea(1000, 5))),
#     type=rep(c("Point", "Area"), nSims)) %>%
#     ggplot(aes(x=prob, group=type, fill=type)) +
#     geom_density(alpha=.2)

# Just get the observed Aproj since thats what we care about for likelihood
AprojObs <- proj$proj$A[simValues$obsField,]

# create matrix for which points fall within areal unit
AprojPoly <- Matrix(
    sparse = T, 
    data = sapply(1:nrow(US.df@data), function(i){
        pix <- rep(0, nrow(AprojObs))
        ones <- simValues %>%
            filter(id==i) %>%
            select(obsAproj) %>%
            unlist %>%
            unname
        pix[ones] <- 1/length(ones)
        pix
    }))

runModel <- function(DFpoint=NULL, DFpoly=NULL, recompile=F, symboic=T, draws=1000){
    model <- "pppSim"
    if(recompile){
        if (file.exists(paste0(model, ".so"))) file.remove(paste0(model, ".so"))
        if (file.exists(paste0(model, ".o"))) file.remove(paste0(model, ".o"))
        if (file.exists(paste0(model, ".dll"))) file.remove(paste0(model, ".dll"))
    }
    compile(paste0(model, ".cpp"))
    if(is.null(DFpoly)){
        empty <- vector("integer") 
        DFpoly <- data.table(obs=empty, denom=empty, loc=empty, long=empty)
    }
    if(is.null(DFpoint)){
        empty <- vector("integer") 
        DFpoint <- data.table(obs=empty, denom=empty, lat=empty, long=empty)
    }
    AprojPoint <- inla.spde.make.A(
        mesh=mesh, 
        loc=as.matrix(select(DFpoint, long, lat)))
    Data <- list(
        yPoint=DFpoint$obs, denomPoint=DFpoint$denom, 
        yPoly=DFpoly$obs, denomPoly=DFpoly$denom, loc=DFpoly$loc,
        M0=spde$param.inla$M0, M1=spde$param.inla$M1, M2=spde$param.inla$M2,
        AprojPoint=AprojPoint, AprojPoly=AprojPoly, AprojObs=AprojObs)
    
    Params <- list(
        beta0=0, log_tau=0, log_kappa=0, z=rep(0, nrow(Q))
    )
    
    dyn.load(dynlib(model))
    startTime <- Sys.time()
    Obj <- MakeADFun(data=Data, parameters=Params, DLL=model, random="z")
    symbolic <- T
    if(symbolic){
        runSymbolicAnalysis(Obj)
    }
    Opt <- nlminb(
        start=Obj$par,
        objective=Obj$fn,
        gradient=Obj$gr,
        control=list(eval.max=1e4, iter.max=1e4))
    sdrep <- sdreport(Obj, getJointPrecision=TRUE)
    runtime <- Sys.time() - startTime
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

N <- 1200
obsPoints <- spsample(US.df, N, "random")@coords
AprojPoint <- inla.spde.make.A(mesh=mesh, loc=obsPoints)
obsDF <- data.table(long=obsPoints[,1], lat=obsPoints[,2]) %>%
    mutate(re=as.vector(AprojPoint %*% x), denom=rpois(N, 100)) %>%
    mutate(prob=invlogit(beta0 + re), obs=rbinom(N, denom, prob))

USDF %>%
    ggplot +
    aes(long,lat,group=group) + 
    geom_path(color="black", size=.1) +
    coord_equal() + 
    theme_void() +
    geom_point(
        aes(group=NULL), 
        data=obsDF,
        color="red",
        size=.5)

plot(mesh)
points(obsPoints, col="red", pch=21, cex=.2)

simPointsOnly <- runModel(obsDF, recompile=F, symboic=F)
estValues <- mesh2DF(simPointsOnly)

rbind(
    mutate(estValues, type="Estimated"), 
    mutate(simValues, type="True")) %>%
    filter(obsField) %>%
    ggplot(aes(x, y, z=obs)) +
    geom_raster(aes(fill = obs)) + 
    coord_equal() +
    theme_void() + 
    scale_fill_distiller(palette = "Spectral") +
    geom_path(
        aes(long,lat, group=group, fill=NULL, z=NULL),
        color="black",
        size=.1,
        data=USDF) +
    facet_wrap(~type)

# now lets try adding in some polygons to the equation
obsPointsMix <- spsample(US.df, N/2, "random")
AprojPointMix <- inla.spde.make.A(mesh=mesh, loc=obsPointsMix@coords)
obsPointMixDF <- data.table(
    long=obsPointsMix@coords[,1], 
    lat=obsPointsMix@coords[,2]) %>%
    mutate(re=as.vector(AprojPointMix %*% x), denom=rpois(n(), 100)) %>%
    mutate(prob=invlogit(beta0 + re), obs=rbinom(n(), denom, prob))

obsPolyMixDF <- over(obsPointsMix, US.df) %>%
    group_by(id) %>%
    summarize(nObs=n()) %>%
    right_join(countyRiskDF, by="id") %>%
    filter(is.na(nObs)) %>%
    select(-nObs) %>%
    rename(re=obs) %>%
    mutate(denom=cells*5, loc=id-1) %>%
    mutate(prob=invlogit(beta0 + re), obs=rbinom(n(), denom, prob))

USDF %>%
    left_join(obsPolyMixDF) %>%
    mutate(observed=as.numeric(!is.na(denom)) * .2) %>%
    ggplot +
    aes(long,lat,group=group) + 
    geom_path(color="black", size=.1) +
    geom_polygon(aes(alpha=observed)) +
    coord_equal() + 
    theme_void() +
    geom_point(
        aes(group=NULL), 
        data=obsPointMixDF,
        color="red",
        size=.5) +
    guides(alpha=FALSE)

simPPMix <- runModel(obsPointMixDF, obsPolyMixDF, recompile=F, symboic=F)
estValuesMix <- mesh2DF(simPPMix)
rbind(
    mutate(estValuesMix, type="Estimated Mix"),
    mutate(estValues, type="Estimated Points"),
    mutate(simValues, type="True")) %>%
    filter(obsField) %>%
    ggplot(aes(x, y, z=obs)) +
    geom_raster(aes(fill = obs)) + 
    coord_equal() +
    theme_void() + 
    scale_fill_distiller(palette = "Spectral") +
    geom_path(
        aes(long,lat, group=group, fill=NULL, z=NULL),
        color="black",
        size=.1,
        data=USDF) +
    facet_wrap(~type, nrow = 2)

## What about only polys

obsPolyDF <- countyRiskDF %>%
    arrange(id) %>%
    rename(re=obs) %>%
    mutate(denom=cells*5, loc=id-1) %>%
    mutate(prob=invlogit(beta0 + re), obs=rbinom(n(), denom, prob))

simPolyOnly <- runModel(NULL, obsPolyDF, recompile=F, symboic=F)
estValuesPoly <- mesh2DF(simPolyOnly)

rbind(
    mutate(estValuesMix, type="Estimated Mix"),
    mutate(estValues, type="Estimated Points"),
    mutate(estValuesPoly, type="Estimated Polygons"),
    mutate(simValues, type="True")) %>%
    filter(obsField) %>%
    ggplot(aes(x, y, z=obs)) +
    geom_raster(aes(fill = obs)) + 
    coord_equal() +
    theme_void() + 
    scale_fill_distiller(palette = "Spectral") +
    geom_path(
        aes(long,lat, group=group, fill=NULL, z=NULL),
        color="black",
        size=.1,
        data=USDF) +
    facet_wrap(~type, nrow = 2)

rbind(
    mutate(estValuesMix, type="Estimated Mix"),
    mutate(estValues, type="Estimated Points"),
    mutate(estValuesPoly, type="Estimated Polygons"),
    mutate(simValues, type="True")) %>%
    filter(obsField) %>%
    group_by(id, type) %>%
    summarize(obs=mean(obs, na.rm=T)) %>%
    full_join(USDF, by="id", copy=T) %>%
    ggplot +
    aes(long,lat,group=group,fill=obs) +
    geom_polygon() +
    geom_path(color="black", size=.1) +
    coord_equal() +
    theme_void() +
    scale_fill_distiller(palette = "Spectral") +
    facet_wrap(~type, nrow = 2)

rbind(
    mutate(estValuesMix, type="Estimated Mix"),
    mutate(estValues, type="Estimated Points"),
    mutate(estValuesPoly, type="Estimated Polygons")) %>%
    filter(obsField) %>%
    ggplot(aes(x, y, z=sdx)) +
    geom_raster(aes(fill = sdx)) + 
    coord_equal() +
    theme_void() + 
    scale_fill_distiller(palette = "Spectral") +
    geom_path(
        aes(long,lat, group=group, fill=NULL, z=NULL),
        color="black",
        size=.1,
        data=USDF) +
    facet_wrap(~type, nrow = 2)

evalModel <- function(model){
    c(rmse=sqrt(mean((model$z$mu - x)^2)),
      coverage=mean((x >= model$z$lwr) & (x <= model$z$upr)),
      runtime=model$runtime)
}

list(Points=simPointsOnly, Poly=simPolyOnly, Mix=simPPMix) %>%
    sapply(evalModel)

rbind(
    mutate(estValuesMix, type="Estimated Mix"),
    mutate(estValues, type="Estimated Points"),
    mutate(estValuesPoly, type="Estimated Polygons")) %>%
    filter(obsField) %>%
    group_by(type) %>%
    summarize(sdxm=mean(sdx))
