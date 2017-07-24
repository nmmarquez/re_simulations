#!/homes/nmarquez/packages/R-3.4.1/bin/Rscript

rm(list=ls())
INLA:::inla.dynload.workaround()
pacman::p_load(INLA, ggplot2, data.table, lattice, TMB, ar.matrix, MASS,
               argparse, clusterPower)
set.seed(123)

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
parser$add_argument("-N", required=TRUE, type="integer",
                    help="Number of spatial points.")

args <- parser$parse_args()

mesh_to_dt <- function(x, proj, time, model){
    M <- length(proj$x)
    DT <- data.table(x=rep(proj$x, M), y=rep(proj$y, each=M), time=time,
                     obs=c(inla.mesh.project(proj, field=x)), model=model)
    DT
}


n <- args$N # number of observations on the grid
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

ParList <- list(rho=rho, kappa=kappa0, tau=tau0, sigma=sigma0, range=range0,
                N=args$N)
save_folder <- "/home/j/temp/nmarquez/sim2dresults/"
outf <- paste0(save_folder, "sigma_", ParList$sigma, "_range_",
               ParList$range, "_rho_", ParList$rho, "_N_", ParList$N,
               "_summary.Rout")
output_file <- file(outf, open="wt")
sink(output_file)
sink(output_file, type="message")

spde <- inla.spde2.matern(mesh) # create the spde from the mesh

# calculate the precision matrix by hand
Qgeo <- tau0**2 * (kappa0**4 * spde$param.inla$M0 +
                   2 * kappa0**2 *spde$param.inla$M1 + spde$param.inla$M2)

# sim by krnoecker
Q <- kronecker(Q.AR1(m, 1, rho), Qgeo)
x_ <- matrix(data=c(sim.AR(1, Q)), nrow=mesh$n, ncol=m)
x_ <- x_ - mean(x_)

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
inla.setOption(num.threads=4) 
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
               phi=array(0, dim=c(nrow(Qgeo),m)))

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

inlares <- sapply(1:m, function(i) res$summary.random$i$mean[iset$i.group==i])
tmbres <- Report$phi
trueres <- x_

(modelabsdiff <- mean(abs(tmbres - inlares)))
(inlabsdiff <- mean(abs(trueres - inlares)))
(tmbabsdiff <- mean(abs(tmbres - trueres)))

Qphi <- sdrep$jointPrecision[row.names(sdrep$jointPrecision) == "phi",
                             row.names(sdrep$jointPrecision) == "phi"]

system.time(phi_draws <- t(sim.AR(1000, Qphi)) + c(tmbres)) # 40 seconds
system.time(draws <- inla.posterior.sample(1000, res)) # 76 seconds
inla_draws <- sapply(1:1000, function(x)
    draws[[x]]$latent[grepl("i:", row.names(draws[[x]]$latent)), 1])
inlabdraws <- sapply(1:1000, function(x) draws[[x]]$latent[c("A", "B", "C"),])

VC <- sdrep$jointPrecision[row.names(sdrep$jointPrecision) != "phi",
                           row.names(sdrep$jointPrecision) != "phi"]
nonbetas <- row.names(VC)[4:nrow(VC)]
tmbfdraws <- t(sim.AR(1000, VC)) + c(Report$beta, sapply(nonbetas, function(x)
    Report[[x]]))
row.names(tmbfdraws) <- row.names(VC)
# exp rows
tvec <- c(beta="beta", logtau="tau", logkappa="kappa",
          logitrho="rho", logsigma="sigma")

for(r in c("logtau", "logkappa", "logsigma")){
    tmbfdraws[r,] <- exp(tmbfdraws[r,])
}

for(r in c("logitrho")){
    tmbfdraws[r,] <- expit(tmbfdraws[r,])
}

row.names(tmbfdraws) <- tvec[row.names(tmbfdraws)]
row.names(tmbfdraws)[1:3] <- paste0("beta", 1:3)

hpdraws <- t(mapply(function(x,y) rnorm(1000, x, y), summary(res)$hyperpar$mean,
                    summary(res)$hyperpar$sd))
hpdraws[1,] <- 1/hpdraws[1,]**.5
hpdraws[2,] <- exp(hpdraws[2,])**.5
hpdraws[3,] <- exp(hpdraws[3,])

DTfixed <- rbindlist(list(
    data.table(value=c(tmbfdraws), draw=rep(1:1000, each=nrow(tmbfdraws)),
               par=rep(row.names(tmbfdraws), 1000), method="TMB"),
    data.table(value=c(inlabdraws), draw=rep(1:1000, each=3),
               par=rep(paste0("beta", 1:3), 1000), method="INLA"),
    data.table(value=c(hpdraws), draw=rep(1:1000, each=4),
               par=rep(c("sigma", "tau", "kappa", "rho"), 1000), method="INLA")))

DTpars <- data.table(value=c(sd.y, tau0, kappa0, rho, -1:1),
                     par=c("sigma", "tau", "kappa", "rho", paste0("beta", 1:3)))


inlavarlist <- lapply(1:m, function(i)
    mesh_to_dt(res$summary.random$i$sd[iset$i.group==i], proj, i, "inla"))
tmbvarlist <- lapply(1:m, function(i)
    mesh_to_dt(apply(phi_draws[iset$i.group==i,],1,sd), proj, i, "tmb"))

DTvar <- rbindlist(c(inlavarlist, tmbvarlist))

tmbpreds <- phi_draws[rep(1:mesh$n %in% mesh$idx$loc, m),] +
    sapply(c("A", "B", "C"), function(a) as.integer(a == ccov)) %*% tmbfdraws[1:3,] +
    t(sapply(1:(n*m), function(x) rnorm(1000, 0, exp(Report$logsigma))))

inlapreds <- inla_draws[rep(1:mesh$n %in% mesh$idx$loc, m),] +
    sapply(c("A", "B", "C"), function(a) as.integer(a == ccov)) %*% inlabdraws +
    t(sapply(1:(n*m), function(x)
        rnorm(1000, 0, (res$summary.hyperpar[1,"mean"])**-.5)))

summary(apply(tmbpreds, 1, mean))
summary(apply(inlapreds, 1, mean))
mean(abs(apply(inlapreds,1, mean) - apply(tmbpreds, 1, mean)))


tmbquant <- t(apply(tmbpreds, 1, quantile, probs=c(.025, .975)))
inlaquant <- t(apply(inlapreds, 1, quantile, probs=c(.025, .975)))

# in sample
(incovtmb <- mean((c(y)[isel] >= tmbquant[isel,1]) & (c(y)[isel] <= tmbquant[isel,2])))
(incovinla <- mean((c(y)[isel] >= inlaquant[isel,1]) & (c(y)[isel] <= inlaquant[isel,2])))
(inrmsetmb <- mean((c(y)[isel] - apply(tmbpreds, 1, mean)[isel])**2)**.5)
(inrmseinla <- mean((c(y)[isel] - apply(inlapreds, 1, mean)[isel])**2)**.5)
(inmadtmb <- median(abs(c(y)[isel] - apply(tmbpreds, 1, mean)[isel])))
(inmadinla <- median(abs(c(y)[isel] - apply(inlapreds, 1, mean)[isel])))
(inrmsediff <- inrmseinla - inrmsetmb)
(inmaddiff <- inmadinla - inmadtmb)


# out of sample
(outcovtmb <- mean((c(y)[-isel] >= tmbquant[-isel,1]) & (c(y)[-isel] <= tmbquant[-isel,2])))
(outcovinla <- mean((c(y)[-isel] >= inlaquant[-isel,1]) & (c(y)[-isel] <= inlaquant[-isel,2])))
(outrmsetmb <- mean((c(y)[-isel] - apply(tmbpreds, 1, mean)[-isel])**2)**.5)
(outrmseinla <- mean((c(y)[-isel] - apply(inlapreds, 1, mean)[-isel])**2)**.5)
(outmadtmb <- median(abs(c(y)[-isel] - apply(tmbpreds, 1, mean)[-isel])))
(outmadinla <- median(abs(c(y)[-isel] - apply(inlapreds, 1, mean)[-isel])))
(outrmsediff <- outrmseinla - outrmsetmb)
(outmaddiff <- outmadinla - outmadtmb)

DataList <- list(variance=DTvar, params=DTpars, fixed=DTfixed, latent=DT)
DataList <- c(ParList, DataList)
MetaList <- list(inrmsediff=inrmsediff, outrmsediff=outrmsediff, 
                 inmaddiff=inmaddiff, outmaddiff=outmaddiff, 
                 inla.time=inla.time, tmb.time=tmb.time, outcovtmb=outcovtmb,
                 incovtmb=incovtmb)
MetaList <- c(ParList, MetaList)

save_file_data <- paste0(save_folder, "sigma_", ParList$sigma, "_range_",
                         ParList$range, "_rho_", ParList$rho, "_N_", ParList$N,
                         "_data.Rda")
save_file_meta <- paste0(save_folder, "sigma_", ParList$sigma, "_range_",
                         ParList$range, "_rho_", ParList$rho, "_N_", ParList$N,
                         "_meta.Rda")
save(DataList, file=save_file_data)
save(MetaList, file=save_file_meta)

sink()
sink(type="message")
