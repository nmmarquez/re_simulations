rm(list=ls())
output_file <- file("~/Documents/re_simulations/inla/sim2d.Rout", open="wt")
sink(output_file)
sink(output_file, type="message")
pacman::p_load(INLA, ggplot2, data.table, lattice, TMB, ar.matrix, MASS)
set.seed(123)
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

n <- 500 # number of observations on the grid
m <- 12 # number of time points
loc <- matrix(runif(n*2), n, 2) # simulate observed points
mesh <- inla.mesh.create(loc, refine=list(max.edge=0.05)) # create mesh
jpeg("~/Documents/re_simulations/inla/mesh.jpg")
par(mfrow=c(1,1))
plot(mesh)
points(loc[,1], loc[,2], col="red", pch=20)
dev.off()
# project mesh using inla default projection
proj <- inla.mesh.projector(mesh)

sigma0 <-  .1   # Standard deviation
range0 <- .5 # Spatial range
kappa0 <- sqrt(8)/range0 # inla paramter transform
tau0 <- 1/(sqrt(4*pi)*kappa0*sigma0) # inla parameter transform
rho <- .91 # tenporal autocorrelation

spde <- inla.spde2.matern(mesh) # create the spde from the mesh

# parameterize the spde and get the projected precision matrix
Q1 <- inla.spde2.precision(spde, theta=c(log(tau0), log(kappa0)))
# calculate the precision matrix by hand in order to make sure you got the
# process down for TMB
Q2 <- tau0**2 * (kappa0**4 * spde$param.inla$M0 + 
                     2 * kappa0**2 *spde$param.inla$M1 + spde$param.inla$M2)
# Should all be equal
print(all.equal(Q1, Q2))

# # simulate m sets from the precision matrix had no idea INLA had this!!!
# x.m <- inla.qsample(n=m, Q1)
# 
# # use the janky sim code found here 
# # http://www.math.ntnu.no/inla/r-inla.org/tutorials/spde/html/
# x_ <- x.m
# for (j in 2:m) 
#     x_[,j] <- rho*x_[,j-1] + sqrt(1-rho^2)*x.m[,j]

# sim by krnoecker
Q <- kronecker(Q.AR1(m, 1, rho), Q1)
x_ <- matrix(data=c(sim.AR(1, Q)), nrow=mesh$n, ncol=m)
x_ <- x_ - mean(x_)

# lets only take the observed mesh not the whole set
x <- x_[mesh$idx$loc,]

# plot using the above websites code
c100 <- rainbow(101)

jpeg("~/Documents/re_simulations/inla/pointplot.jpg")
par(mfrow=c(4,3), mar=c(0,0,0,0))
for (j in 1:m)
    plot(loc, col=c100[round(100*(x[,j]-min(x[,j]))/diff(range(x[,j])))], 
         axes=FALSE, asp=1, pch=19, cex=0.5)
dev.off()

# plot using our use defined code to see the whole surface
for (j in 1:m){
   print(plot_mesh_sim(x_[,j], proj) + labs(title=paste0("Time: ", j)))
}

# lets build up a linear model with dummies to estimate
table(ccov <- factor(sample(LETTERS[1:3], n*m, replace=TRUE)))

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
print(system.time(res <- inla(formulae,  data=inla.stack.data(sdat),
                    control.predictor=list(compute=TRUE, A=inla.stack.A(sdat)),
                    control.family=list(hyper=list(theta=prec.prior)),
                    control.fixed=list(expand.factor.strategy='inla'),
                    control.compute=list(config = TRUE),
                    control.inla=list(int.strategy="eb"))))
# 2440(7018.576) seconds run time
summary(res)


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
Obj <- MakeADFun(data=Data, parameters=Params, DLL=model, random="phi")
runSymbolicAnalysis(Obj)
print(system.time(Opt <- nlminb(start=Obj$par, objective=Obj$fn,
                                gradient=Obj$gr,
                                control=list(eval.max=1e4, iter.max=1e4))))
# get the estimated values
Report <- Obj$report()
system.time(sdrep <- sdreport(Obj, getJointPrecision = T))
Q <- sdrep$jointPrecision[row.names(sdrep$jointPrecision) == "phi", 
                          row.names(sdrep$jointPrecision) == "phi"]

# 573.485(858.640) + 89.754(247.460) = 663.239(1106.100)

# save the results
save(list=ls(), file="~/Documents/re_simulations/inla/model_results.Rda")

sink(type="message")
sink()