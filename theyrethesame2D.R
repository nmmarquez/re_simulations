rm(list=ls())
pacman::p_load(INLA, ggplot2, data.table, lattice, TMB)
set.seed(123)
# compare the INLA Q matrix vs the by hand to make sure we are on the same page

plot_mesh_sim <- function(x, proj){
    M <- length(proj$x)
    DT <- data.table(x=rep(proj$x, M), y=rep(proj$y, each=M), 
                     obs=c(inla.mesh.project(proj, field=x)))
    ggplot(DT, aes(x, y, z= obs)) + geom_tile(aes(fill = obs)) + theme_bw() + 
        lims(y=c(0,1), x=c(0,1)) + scale_fill_gradientn(colors=heat.colors(8))
}

n <- 500
m <- 12
loc <- matrix(runif(n*2), n, 2)
mesh <- inla.mesh.create(loc, refine=list(max.edge=0.05))
par(mfrow=c(1,1))
plot(mesh)
points(loc[,1], loc[,2], col="red", pch=20)
proj <- inla.mesh.projector(mesh)

sigma0 <-  .3   ## Standard deviation
range0 <- 1. ## Spatial range
kappa0 <- sqrt(8)/range0
tau0 <- 1/(sqrt(4*pi)*kappa0*sigma0)
rho <- .91

spde <- inla.spde2.matern(mesh)#, prior.range=c(0.5, 0.01), prior.sigma=c(1, 0.01))
Q1 <- inla.spde2.precision(spde, theta=c(log(tau0), log(kappa0)))
Q2 <- tau0**2 * (kappa0**4 * spde$param.inla$M0 + 
                     2 * kappa0**2 *spde$param.inla$M1 + spde$param.inla$M2)
all.equal(Q1, Q2)
x.m <- inla.qsample(n=m, Q1)

x_ <- x.m
for (j in 2:m) 
    x_[,j] <- rho*x_[,j-1] + sqrt(1-rho^2)*x.m[,j]

x <- x_[mesh$idx$loc,]

c100 <- rainbow(101)
par(mfrow=c(4,3), mar=c(0,0,0,0))
for (j in 1:m)
    plot(loc, col=c100[round(100*(x[,j]-min(x[,j]))/diff(range(x[,j])))], 
         axes=FALSE, asp=1, pch=19, cex=0.5)

for (j in 1:m){
    print(plot_mesh_sim(x_[,j], proj))
}

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
                  geo=mesh$idx$loc)[isel, ] 

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
# res <- inla(formulae,  data=inla.stack.data(sdat), 
#             control.predictor=list(compute=TRUE, A=inla.stack.A(sdat)), 
#             control.family=list(hyper=list(theta=prec.prior)), 
#             control.fixed=list(expand.factor.strategy='inla'))
# 
# summary(res)

setwd("~/Documents/re_simulations/inla/")

model <- "st"
compile(paste0(model, ".cpp"))

Data <- list(y=y, T=m, geo=dat$geo, temporal=dat$time,
             cov=sapply(c("A", "B", "C"), function(x) as.integer(x == dat$w)),
             M0=spde$param.inla$M0, M1=spde$param.inla$M1, M2=spde$param.inla$M2)
Params <- list(logtau=0, logsigma=0, logitrho=0, logkappa=0, beta=c(0,0,0),
               phi=array(0, dim=c(nrow(Q1),m)))

dyn.load(dynlib(model))
Obj <- MakeADFun(data=Data, parameters=Params, DLL=model, random="phi",
                 silent=!verbose)
Obj$env$tracemgc <- verbose
Obj$env$inner.control$trace <- verbose
print(system.time(Opt <- nlminb(start=Obj$par, objective=Obj$fn, 
                                gradient=Obj$gr,
                                control=list(eval.max=1e4, iter.max=1e4))))
Report <- Obj$report()