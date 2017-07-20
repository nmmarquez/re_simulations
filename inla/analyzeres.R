rm(list=ls())
pacman::p_load(data.table, ggplot2, INLA, TMB, ar.matrix, clusterPower)
load(file="~/Documents/re_simulations/inla/model_results.Rda")

mesh_to_dt <- function(x, proj, time, model){
    M <- length(proj$x)
    DT <- data.table(x=rep(proj$x, M), y=rep(proj$y, each=M), time=time, 
                     obs=c(inla.mesh.project(proj, field=x)), model=model)
    DT
}

datalist <- lapply(1:m, function(i) 
    mesh_to_dt(x_[,i] - mean(x_), proj, i, "data"))
inlalist <- lapply(1:m, function(i) 
    mesh_to_dt(res$summary.random$i$mean[iset$i.group==i], proj, i, "inla"))
tmblist <- lapply(1:m, function(i) 
    mesh_to_dt(Report$phi[,i], proj, i, "tmb"))


DT <- rbindlist(c(datalist, inlalist, tmblist))

# for(i in 1:m){
#     print(ggplot(DT[time==i,], aes(x, y, z= obs)) + geom_tile(aes(fill = obs)) + 
#           theme_bw() + lims(y=c(0,1), x=c(0,1)) + facet_wrap(~model) +
#           scale_fill_gradientn(colors=heat.colors(8)) +
#           labs(title=paste0("Time Point: ", i)))
# }

ggplot(DT[time %in% 1:3,], aes(x, y, z= obs)) + geom_tile(aes(fill = obs)) + 
    theme_bw() + lims(y=c(0,1), x=c(0,1)) + facet_grid(model~time) +
    scale_fill_gradientn(colors=heat.colors(8))

c(sd.y, 1 / res$summary.hyperpar[1,"mean"]**.5, exp(Report$logsigma))
c(tau0, exp(res$summary.hyperpar[2,"mean"])**.5, exp(Report$logtau))
c(kappa0, exp(c(res$summary.hyperpar[3,"mean"], Report$logkappa)))
c(rho, res$summary.hyperpar[4,"mean"], Report$rho)

inlares <- sapply(1:m, function(i) res$summary.random$i$mean[iset$i.group==i])
tmbres <- Report$phi
trueres <- sapply(1:m, function(i) x_[,i] - mean(x_[,i]))

mean(abs(tmbres - inlares))
mean(abs(trueres - inlares))
mean(abs(tmbres - trueres))

max(abs(tmbres - inlares))
max(abs(trueres - inlares))
max(abs(tmbres - trueres))

Qphi <- sdrep$jointPrecision[row.names(sdrep$jointPrecision) == "phi", 
                             row.names(sdrep$jointPrecision) == "phi"]

#system.time(phi_draws <- inla.qsample(1000, Qphi) + c(tmbres)) # 70 seconds
system.time(phi_draws <- t(sim.AR(1000, Qphi)) + c(tmbres)) # 40 seconds
system.time(draws <- inla.posterior.sample(1000, res)) # 400 seconds
inla_draws <- sapply(1:1000, function(x) 
    draws[[x]]$latent[grepl("i:", row.names(draws[[x]]$latent)), 1])
inlabdraws <- sapply(1:1000, function(x) draws[[x]]$latent[c("A", "B", "C"),])
apply(inlabdraws, 1, mean)

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


ggplot(data=DTfixed, aes(x=value, fill=method, group=method)) + 
    geom_density(alpha=.6) + facet_wrap(~par, scales="free")
ggplot(data=DTfixed[par=="kappa"], aes(x=value, fill=method, group=method)) + 
    geom_density()
ggplot(data=DTfixed[par=="beta1"], aes(x=value, fill=method, group=method)) + 
    geom_density()
apply(tmbfdraws, 1, mean)
apply(tmbfdraws, 1, sd)



inla_bounds <- t(apply(inla_draws, 1, quantile, probs=c(.025, .975)))
tmb_bounds <- t(apply(phi_draws, 1, quantile, probs=c(.025, .975)))

mean(c(trueres >= inla_bounds[,1] & trueres <= inla_bounds[,2]))
mean(c(trueres >= tmb_bounds[,1] & trueres <= tmb_bounds[,2]))


inla_sd_draw <- matrix(apply(inla_draws, 1, sd), nrow=mesh$n ,ncol=m)
inla_sd_direct <- sapply(1:m, function(i) 
    res$summary.random$i$sd[iset$i.group==i])
mean(abs(inla_sd_direct - inla_sd_draw))



inlavarlist <- lapply(1:m, function(i)
    mesh_to_dt(res$summary.random$i$sd[iset$i.group==i], proj, i, "inla"))
tmbvarlist <- lapply(1:m, function(i)
    mesh_to_dt(apply(phi_draws[iset$i.group==i,],1,sd), proj, i, "tmb"))

DTvar <- rbindlist(c(inlavarlist, tmbvarlist))
ggplot(DTvar[obs <=.72 & model=="inla"], aes(x, y, z= obs)) + geom_tile(aes(fill = obs)) + 
    theme_bw() + lims(y=c(-.15,1.15), x=c(-.15,1.15)) + facet_wrap(~time) +
    scale_fill_gradientn(colors=heat.colors(8))
ggplot(DTvar[obs <=.15 & model=="tmb"], aes(x, y, z= obs)) + geom_tile(aes(fill = obs)) + 
    theme_bw() + lims(y=c(-.15,1.15), x=c(-.15,1.15)) + facet_wrap(~time) +
    scale_fill_gradientn(colors=heat.colors(8))

ggplot(DTvar[obs & model=="inla"], aes(x, y, z= obs)) + geom_tile(aes(fill = obs)) + 
    theme_bw() + lims(y=c(-.15,1.15), x=c(-.15,1.15)) + facet_wrap(~time) +
    scale_fill_gradientn(colors=heat.colors(8))
ggplot(DTvar[obs & model=="tmb"], aes(x, y, z= obs)) + geom_tile(aes(fill = obs)) + 
    theme_bw() + lims(y=c(-.15,1.15), x=c(-.15,1.15)) + facet_wrap(~time) +
    scale_fill_gradientn(colors=heat.colors(8))
