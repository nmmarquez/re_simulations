rm(list=ls())
pacman::p_load(data.table, ggplot2, INLA, TMB, ar.matrix)
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

c(sd.y, 1 / res$summary.hyperpar[1,"mean"]**.5, Report$sigma)
c(tau0, exp(c(res$summary.hyperpar[2,"mean"], Report$logtau))**.5)
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
