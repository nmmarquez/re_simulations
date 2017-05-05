rm(list=ls())
pacman::p_load(INLA, ggplot2, data.table, lattice)
set.seed(123)
# compare the INLA Q matrix vs the by hand to make sure we are on the same page

plot_mesh_sim <- function(x, proj){
    M <- length(proj$x)
    DT <- data.table(x=rep(proj$x, M), y=rep(proj$y, each=M), 
                     obs=c(inla.mesh.project(proj, field=x)))
    ggplot(DT, aes(x, y, z= obs)) + geom_tile(aes(fill = obs)) + theme_bw() + 
        lims(y=c(0,1), x=c(0,1)) + scale_fill_gradientn(colors=heat.colors(8))
}

n <- 1000
loc <- matrix(runif(n*2), n, 2)
mesh <- inla.mesh.create(loc, refine=list(max.edge=0.05))
plot(mesh)
points(loc[,1], loc[,2], col="red", pch=20)
proj <- inla.mesh.projector(mesh)

sigma0 <-  .3   ## Standard deviation
range0 <- 1. ## Spatial range
kappa0 <- sqrt(8)/range0
tau0 <- 1/(sqrt(4*pi)*kappa0*sigma0)

spde <- inla.spde2.matern(mesh)

##############
Q1 <- inla.spde2.precision(spde, theta=c(log(tau0), log(kappa0)))
Q2 <- tau0**2 * (kappa0**4 * spde$param.inla$M0 + 
                     2 * kappa0**2 *spde$param.inla$M1 + spde$param.inla$M2)
x <- as.vector(inla.qsample(n=1, Q2))

plot_mesh_sim(x, proj)
cov1 <- rnorm(n)

y <- 3 + 2 * cov1 + x[mesh$idx$loc] + rnorm(n, sd=.1)

DT <- data.table(y=y, cov=cov1, id=1:n)

spatform <- y ~ cov1 + f(id, model=spde)

system.time(spatmodel <- inla(spatform, data=DT))

summary(spatmodel)
spatialhat <- inla.spde2.result(spatmodel, "id", spde)
xhat <- spatialhat$summary.values["0.5quant"]
plot_mesh_sim(c(xhat)[[1]], proj)

exp(spatmodel$summary.hyperpar$mean[2:3])
c(tau0, kappa0)
