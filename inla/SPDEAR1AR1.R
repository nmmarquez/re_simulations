rm(list=ls())
output_file <- file("~/Documents/re_simulations/inla/sim3d.Rout", open="wt")
sink(output_file)
sink(output_file, type="message")
pacman::p_load(INLA, ggplot2, data.table, lattice, TMB, ar.matrix)
set.seed(123)

# We are in 3d territory now without INLA so lets jut run the TMB model

# plotting utility function
plot_mesh_sim <- function(x, proj){
    M <- length(proj$x)
    DT <- data.table(x=rep(proj$x, M), y=rep(proj$y, each=M), 
                     obs=c(inla.mesh.project(proj, field=x)))
    ggplot(DT, aes(x, y, z= obs)) + geom_tile(aes(fill = obs)) + theme_bw() + 
        lims(y=c(0,1), x=c(0,1)) + scale_fill_gradientn(colors=heat.colors(8))
}

n <- 100 # number of observations on the grid
m <- 12 # number of time points
k <- 6 # number of ages
loc <- matrix(runif(n*2), n, 2) # simulate observed points
mesh <- inla.mesh.create(loc, refine=list(max.edge=1.)) # create mesh
jpeg("~/Documents/re_simulations/inla/mesh3d.jpg")
par(mfrow=c(1,1))
plot(mesh)
points(loc[,1], loc[,2], col="red", pch=20)
dev.off()
# project mesh using inla default projection
proj <- inla.mesh.projector(mesh)

sigma0 <-  .3   # Standard deviation
range0 <- 1. # Spatial range
kappa0 <- sqrt(8)/range0 # inla paramter transform
tau0 <- 1/(sqrt(4*pi)*kappa0*sigma0) # inla parameter transform
rho <- c(.96, .85) # temporal and age correlation autocorrelation

spde <- inla.spde2.matern(mesh) # create the spde from the mesh

# calculate the precision matrix by hand in order to make sure you got the
# process down for TMB
Qspde <- tau0**2 * (kappa0**4 * spde$param.inla$M0 + 
                        2 * kappa0**2 *spde$param.inla$M1 + spde$param.inla$M2)

# sim by krnoecker
Q <- kronecker(kronecker(Q.AR1(m, 1, rho[1]), Q.AR1(k, 1, rho[2])), Qspde)
system.time(x_ <- array(data=c(sim.AR(1, Q)), dim=c(mesh$n, k, m)))

# lets only take the observed mesh not the whole set
x <- x_[mesh$idx$loc,,]

# lets build up a linear model with dummies to estimate
table(ccov <- factor(sample(LETTERS[1:3], n*m*k, replace=TRUE)) )

beta <- -1:1

sd.y <- 0.1
y <- beta[unclass(ccov)] + c(x) + rnorm(n*m*k, 0, sd.y)
tapply(y, ccov, mean)

DT <- data.table(expand.grid(geo=mesh$idx$loc-1, age=0:(k-1), time=0:(m-1)))
DT[,w:=ccov]
DT[,y:=as.vector(y)]

isel <- sample(1:(n*m*k), n*m*k/1.3)
DT <- DT[isel, ]

plot_mesh_sim(x_[,1,1], proj)

# now run the TMB model using ./st.cpp
setwd("~/Documents/re_simulations/inla/")

# compile the code if not there
model <- "sta"
#if (file.exists(paste0(model, ".so"))) file.remove(paste0(model, ".so"))
#if (file.exists(paste0(model, ".o"))) file.remove(paste0(model, ".o"))
#if (file.exists(paste0(model, ".dll"))) file.remove(paste0(model, ".dll"))
compile(paste0(model, ".cpp"))

# set the data
Data <- list(y=DT$y, T=m, geo=DT$geo, temporal=DT$time, A=k, age=DT$age,
             cov=sapply(c("A", "B", "C"), function(x) as.integer(x == DT$w)),
             M0=spde$param.inla$M0, M1=spde$param.inla$M1, M2=spde$param.inla$M2)

# set the param starting points
Params <- list(logtau=0, logsigma=0, logitrho=c(0, 0), logkappa=0, 
               beta=c(0,0,0), phi=array(0, dim=c(mesh$n, k, m)))

# load and optimize
dyn.load(dynlib(model))
Obj <- MakeADFun(data=Data, parameters=Params, DLL=model, random="phi")
print(system.time(Opt <- nlminb(start=Obj$par, objective=Obj$fn,
                                gradient=Obj$gr,
                                control=list(eval.max=1e4, iter.max=1e4))))
# get the estimated values
Report <- Obj$report()

# save the results
save(list=ls(), file="~/Documents/re_simulations/inla/model_results3D.Rda")

sink(type="message")
sink()