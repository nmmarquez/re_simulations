rm(list=ls())
pacman::p_load(INLA, ggplot2, data.table, lattice, TMB, ar.matrix, MASS,
               argparse)

# create parser object
parser <- ArgumentParser()
# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("--sigma", required=TRUE, type="double",
                    help="Error from spatial process.")
parser$add_argument("--range", required=TRUE, type="double",
                    help="Spatial Range")
parser$add_argument("--rho", required=TRUE, type="double",
                    help="Temporal Auto-correlation")
parser$add_argument("--seed", required=TRUE, type="integer",
                    help="The random seed state of the R process.")

args <- parser$parse_args()

set.seed(args$seed)
# compare the INLA Q matrix vs the by hand to make sure we are on the same page
# use inla to create the model and make sure results match TMB

# plotting utility function
plot_mesh_sim <- function(x, proj){
    M <- length(proj$x)
    DT <- data.table(x=rep(proj$x, M), y=rep(proj$y, each=M), 
                     obs=c(inla.mesh.project(proj, field=x)))
    ggplot(DT, aes(x, y, z= obs)) + geom_tile(aes(fill = obs)) + theme_bw() + 
        lims(y=c(0,1), x=c(0,1)) + scale_fill_gradientn(colors=heat.colors(8))
}

mesh_to_dt <- function(x, proj, time, model){
    M <- length(proj$x)
    DT <- data.table(x=rep(proj$x, M), y=rep(proj$y, each=M), time=time, 
                     obs=c(inla.mesh.project(proj, field=x)), model=model)
    DT
}

n <- 500 # number of observations on the grid
m <- 12 # number of time points
loc <- matrix(runif(n*2), n, 2) # simulate observed points
mesh <- inla.mesh.create(loc, refine=list(max.edge=0.05)) # create mesh
# project mesh using inla default projection
proj <- inla.mesh.projector(mesh)

sigma0 <-  round(args$sigma, 2)   # Standard deviation
range0 <- round(args$range, 2) # Spatial range
kappa0 <- sqrt(8)/range0 # inla paramter transform
tau0 <- 1/(sqrt(4*pi)*kappa0*sigma0) # inla parameter transform
rho <- round(args$rho, 2) # tenporal autocorrelation

spde <- inla.spde2.matern(mesh) # create the spde from the mesh

# calculate the precision matrix by hand
Qgeo <- tau0**2 * (kappa0**4 * spde$param.inla$M0 + 
                   2 * kappa0**2 *spde$param.inla$M1 + spde$param.inla$M2)

# sim by krnoecker
Q <- kronecker(Q.AR1(m, 1, rho), Qgeo)
x_ <- matrix(data=c(sim.AR(1, Q)), nrow=mesh$n, ncol=m)

# lets only take the observed mesh not the whole set
x <- x_[mesh$idx$loc,]

# lets build up a linear model with dummies to estimate
table(ccov <- factor(sample(LETTERS[1:3], n*m, replace=TRUE)) )

beta <- -1:1

sd.y <- 0.1
y <- beta[unclass(ccov)] + x + rnorm(n*m, 0, sd.y)
tapply(y, ccov, mean)

isel <- sample(1:(n*m), n*m/2)
dat <- data.table(y=as.vector(y), w=ccov, 
                  time=rep(1:m, each=n) - 1, 
                  xcoo=rep(loc[,1], m), 
                  ycoo=rep(loc[,2], m),
                  geo=mesh$idx$loc-1)[isel, ] 

# build the inla stuff that Im not 100 % sure what it does but theres a 
# projector matrix and some priors and the most stacked lists uve ever seen
iset <- inla.spde.make.index('i', n.spde=spde$n.spde, n.group=m)
A <- inla.spde.make.A(mesh=mesh, 
                      loc=cbind(dat$xcoo, dat$ycoo), 
                      group=dat$time+1) 

sdat <- inla.stack(tag='stdata', data=list(y=dat$y), 
                   A=list(A,1),  effects=list(iset, w=dat$w)) 
h.spec <- list(theta=list(prior='pccor1', param=c(0, 0.9)))
formulae <- y ~ 0 + w + 
    f(i, model=spde, group=i.group, 
      control.group=list(model='ar1', hyper=h.spec)) 
prec.prior <- list(prior='pc.prec', param=c(1, 0.01))

# Run the inla model and time it
start.time <- Sys.time()
res <- inla(formulae,  data=inla.stack.data(sdat),
            control.predictor=list(compute=TRUE, A=inla.stack.A(sdat)),
            control.family=list(hyper=list(theta=prec.prior)),
            control.fixed=list(expand.factor.strategy='inla'),
            control.compute=list(config = TRUE),
            control.inla=list(int.strategy="eb"))
end.time <- Sys.time()
inla.time <- end.time - start.time

# now run the TMB model using ./st.cpp
setwd("~/Documents/re_simulations/inla/")

# compile the code if not there
model <- "st"
if (file.exists(paste0(model, ".so"))) file.remove(paste0(model, ".so"))
if (file.exists(paste0(model, ".o"))) file.remove(paste0(model, ".o"))
if (file.exists(paste0(model, ".dll"))) file.remove(paste0(model, ".dll"))
compile(paste0(model, ".cpp"))

# set the data
Data <- list(y=dat$y, T=m, geo=dat$geo, temporal=dat$time,
             cov=sapply(c("A", "B", "C"), function(x) as.integer(x == dat$w)),
             M0=spde$param.inla$M0, M1=spde$param.inla$M1, M2=spde$param.inla$M2)

# set the param starting points
Params <- list(logtau=0, logsigma=0, logitrho=0, logkappa=0, beta=c(0,0,0),
               phi=array(0, dim=c(nrow(Q1),m)))

# load and optimize
dyn.load(dynlib(model))
start.time <- Sys.time()
Obj <- MakeADFun(data=Data, parameters=Params, DLL=model, random="phi")
runSymbolicAnalysis(Obj)
print(system.time(Opt <- nlminb(start=Obj$par, objective=Obj$fn,
                                gradient=Obj$gr,
                                control=list(eval.max=1e4, iter.max=1e4))))
# get the estimated values
Report <- Obj$report()
system.time(sdrep <- sdreport(Obj, getJointPrecision = T))
Qest <- sdrep$jointPrecision[row.names(sdrep$jointPrecision) == "phi", 
                             row.names(sdrep$jointPrecision) == "phi"]
end.time <- Sys.time()
tmb.time <- end.time - start.time

datalist <- lapply(1:m, function(i) 
    mesh_to_dt(x_[,i] - mean(x_), proj, i, "data"))
inlalist <- lapply(1:m, function(i) 
    mesh_to_dt(res$summary.random$i$mean[iset$i.group==i], proj, i, "inla"))
tmblist <- lapply(1:m, function(i) 
    mesh_to_dt(Report$phi[,i], proj, i, "tmb"))


DT <- rbindlist(c(datalist, inlalist, tmblist))
