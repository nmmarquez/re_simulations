rm(list=ls())
pacman::p_load(INLA, ggplot2, data.table, lattice, arm, dplyr, TMB, ar.matrix)
set.seed(124)

mesh2DF <- function(x){
    M <- length(proj$x)
    DT <- data.table(x=rep(proj$x, M), y=rep(proj$y, each=M), 
                     obs=c(inla.mesh.project(proj, field=x)))
    DT
}

# what is the shape that we are dealing with?
US.df$id <- row.names(US.df@data)
USDF <- fortify(US.df, region="id") %>%
    left_join(US.df@data, by="id")
plot(US.df)

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

simValues <- mesh2DF(x)

simValues %>%
    select(x, y) %>%
    as.matrix %>%
    SpatialPoints(US.df@proj4string) %>%
    over(US.df) %>%
    cbind(simValues) %>%
    filter(!is.na(NAME)) %>%
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
simValues %>%
    select(x, y) %>%
    as.matrix %>%
    SpatialPoints(US.df@proj4string) %>%
    over(US.df) %>%
    cbind(simValues) %>%
    filter(!is.na(NAME)) %>%
    group_by(id) %>%
    summarize(obs=mean(obs, na.rm=T)) %>%
    right_join(USDF) %>%
    ggplot +
    aes(long,lat,group=group,fill=obs) + 
    geom_polygon() +
    geom_path(color="black", size=.1) +
    coord_equal() + 
    theme_void() +
    scale_fill_distiller(palette = "Spectral")
