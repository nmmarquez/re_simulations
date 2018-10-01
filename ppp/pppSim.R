rm(list=ls())
pacman::p_load(INLA, ggplot2, data.table, lattice, arm, dplyr, TMB, 
               ar.matrix, dtplyr)
set.seed(124)

mesh2DF <- function(x){
    M <- length(proj$x)
    DT <- data.table(x=rep(proj$x, M), y=rep(proj$y, each=M), 
                     obs=c(inla.mesh.project(proj, field=x)))
    DT
}

# what is the shape that we are dealing with?
row.names(US.df@data) <- 1:nrow(US.df)
US.df$id <- 1:nrow(US.df)#row.names(US.df@data)
USDF <- fortify(US.df, region="id") %>%
    mutate(id=as.numeric(id)) %>%
    left_join(US.df@data, by="id")

randomSPDF <- spsample(US.df, 800, "random")

USDF %>%
    ggplot +
    aes(long,lat,group=group) + 
    geom_path(color="black", size=.1) +
    coord_equal() + 
    theme_void() +
    geom_point(
        aes(x, y, group=NULL), 
        data=as.data.table(randomSPDF@coords),
        color="red",
        size=.5)

west <- c(
    "El Paso", "Hudspeth", "Culberson", "Reeves", "Pecos", "Terrell", 
    "Jeff Davis", "Presidio", "Brewster")

# check which points fall in west coast

randomSPDF$west <- over(randomSPDF, US.df)$NAME %in% west
randomSPDF$long <- randomSPDF@coords[,1]
randomSPDF$lat <- randomSPDF@coords[,2]

USDF %>%
    filter(NAME %in% west) %>%
    ggplot +
    aes(long,lat,group=group) + 
    geom_path(color="black", size=.1) +
    coord_equal() + 
    theme_void() +
    geom_point(
        aes(long, lat, group=NULL), 
        data=subset(randomSPDF@data, west),
        color="red",
        size=.5)

plot(mesh <- inla.mesh.2d(
    randomSPDF, 
    cutoff=.3,
    max.edge=c(50, 500)))
points(randomSPDF, col="red", pch=21, cex=.3)
proj <- inla.mesh.projector(mesh, dims=c(400, 400))

sigma0 <-  .2   ## Standard deviation
range0 <- 2.5 ## Spatial range
kappa0 <- sqrt(8) / range0
tau0 <- 1/(sqrt(4*pi)*kappa0*sigma0)
spde <- inla.spde2.matern(mesh)

Q <- tau0**2 * (kappa0**4 * spde$param.inla$M0 + 
                     2 * kappa0**2 *spde$param.inla$M1 + spde$param.inla$M2)
x_ <- as.vector(sim.AR(n=1, Q))
x <- x_ - mean(x_)

simValues <- expand.grid(proj$x, proj$y) %>%
    SpatialPoints(US.df@proj4string) %>%
    over(US.df) %>%
    as.data.table %>%
    cbind(mesh2DF(x)) %>%
    mutate(aproj=1:n(), obsField=!is.na(id)) %>%
    group_by(obsField) %>%
    mutate(obsAproj=1:n()) %>%
    ungroup %>%
    as_tibble %>%
    mutate(stateID=group_indices(., STATEFP)) %>%
    mutate(stateID=ifelse(!obsField, NA, stateID)) %>%
    mutate(obsAproj=ifelse(!obsField, NA, obsAproj)) %>%
    as.data.table

simValues %>%
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
        data=USDF)

# assume equal pop weight across pixels
countyRiskDF <- simValues %>%
    filter(obsField) %>%
    group_by(id) %>%
    summarize(obs=mean(obs, na.rm=T)) %>%
    right_join(USDF, by="id", copy=T)

countyRiskDF %>%
    ggplot +
    aes(long,lat,group=group,fill=obs) + 
    geom_polygon() +
    geom_path(color="black", size=.1) +
    coord_equal() + 
    theme_void() +
    scale_fill_distiller(palette = "Spectral")

# Jsut get the observed Aproj since thats what we care about for likelihood
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

model <- "pppSim"
recompile <- TRUE
if(recompile){
    if (file.exists(paste0(model, ".so"))) file.remove(paste0(model, ".so"))
    if (file.exists(paste0(model, ".o"))) file.remove(paste0(model, ".o"))
    if (file.exists(paste0(model, ".dll"))) file.remove(paste0(model, ".dll"))
}

compile(paste0(model, ".cpp"))

N <- 800
beta0 <- -1
obsPoints <- spsample(US.df, N, "random")@coords

AprojPoint <- inla.spde.make.A(mesh=mesh, loc=obsPoints)

obsDF <- data.table(long=obsPoints[,1], lat=obsPoints[,2]) %>%
    mutate(re=as.vector(AprojPoint %*% x), denom=rpois(N, 100)) %>%
    mutate(prob=invlogit(beta0 + re), obs=rbinom(N, denom, prob))

Data <- list(
    yPoint=obsDF$obs, denomPoint=obsDF$denom, 
    yPoly=rep(0,0), denomPoly=rep(0,0),
    M0=spde$param.inla$M0, M1=spde$param.inla$M1, M2=spde$param.inla$M2,
    AprojPoint=AprojPoint, AprojPoly=AprojPoly, AprojObs=AprojObs)

Params <- list(
    beta0=0, log_tau=0, log_kappa=0, z=rep(0, nrow(Q))
)

dyn.load(dynlib(model))
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
Report <- Obj$report()

estValues <- expand.grid(proj$x, proj$y) %>%
    SpatialPoints(US.df@proj4string) %>%
    over(US.df) %>%
    as.data.table %>%
    cbind(mesh2DF(Report$z)) %>%
    mutate(aproj=1:n(), obsField=!is.na(id)) %>%
    group_by(obsField) %>%
    mutate(obsAproj=1:n()) %>%
    ungroup %>%
    as_tibble %>%
    mutate(stateID=group_indices(., STATEFP)) %>%
    mutate(stateID=ifelse(!obsField, NA, stateID)) %>%
    mutate(obsAproj=ifelse(!obsField, NA, obsAproj)) %>%
    as.data.table

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
