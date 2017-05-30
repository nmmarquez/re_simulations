rm(list=ls())
pacman::p_load(data.table, ggplot2, INLA, TMB)
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

ggplot(DT[model=="data",], aes(x, y, z= obs)) + geom_tile(aes(fill = obs)) + 
          theme_bw() + lims(y=c(0,1), x=c(0,1)) + facet_wrap(~model) +
          scale_fill_gradientn(colors=heat.colors(8))

for(i in 1:m){
    print(ggplot(DT[time==i,], aes(x, y, z= obs)) + geom_tile(aes(fill = obs)) + 
          theme_bw() + lims(y=c(0,1), x=c(0,1)) + facet_wrap(~model) +
          scale_fill_gradientn(colors=heat.colors(8)) +
          labs(title=paste0("Time Point: ", i)))
}

ggplot(DT[time %in% 1:3,], aes(x, y, z= obs)) + geom_tile(aes(fill = obs)) + 
    theme_bw() + lims(y=c(0,1), x=c(0,1)) + facet_grid(model~time) +
    scale_fill_gradientn(colors=heat.colors(8))

c(sd.y, 1 / res$summary.hyperpar[1,"mean"]**.5, Report$sigma)
c(tau0, exp(c(res$summary.hyperpar[2,"mean"], Report$logtau))**.5)
c(kappa0, exp(c(res$summary.hyperpar[3,"mean"], Report$logkappa)))
c(rho, res$summary.hyperpar[4,"mean"], Report$rho)

inlares <- sapply(1:m, function(i) res$summary.random$i$mean[iset$i.group==i])
tmbres <- Report$phi
res <- sapply(1:m, function(i) x_[,i] - mean(x_[,i]))

mean(abs(tmbres - inlares))
mean(abs(res - inlares))
mean(abs(tmbres - res))

max(abs(tmbres - inlares))
max(abs(res - inlares))
max(abs(tmbres - res))
