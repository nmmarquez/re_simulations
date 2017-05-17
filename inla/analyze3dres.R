rm(list=ls())
pacman::p_load(data.table, ggplot2, INLA, TMB, ggthemes)
load(file="~/Documents/re_simulations/inla/model_results3D.Rda")

mesh_to_dt <- function(x, proj, age, time, model){
    M <- length(proj$x)
    DT <- data.table(x=rep(proj$x, M), y=rep(proj$y, each=M), time=time+2000, 
                     age=paste(age, "year old"),
                     obs=c(inla.mesh.project(proj, field=x)), model=model)
    DT
}

tdems <- expand.grid(age=1:k, time=1:m)
x_mean <- mean(x_)

datalist <- lapply(1:nrow(tdems), function(i) 
    mesh_to_dt(x_[,tdems$age[i],tdems$time[i]] - x_mean, proj, 
               tdems$age[i], tdems$time[i], "data"))
tmblist <- lapply(1:nrow(tdems), function(i) 
    mesh_to_dt(Report$phi[,tdems$age[i],tdems$time[i]], proj, 
               tdems$age[i], tdems$time[i], "tmb"))


DT <- rbindlist(c(datalist, tmblist))

c(sd.y, Report$sigma)
c(tau0, exp(Report$logtau))
c(kappa0, exp(Report$logkappa))
c(rho, Report$rho)
    
ys_ <- 2001:2004
as_ <- paste(1:2, "year old")

jpeg("~/Documents/re_simulations/inla/proj3drez.jpg")
print(ggplot(DT[time %in% ys_ & age %in% as_,], aes(x, y, z= obs)) +
          geom_tile(aes(fill = obs)) + lims(y=c(0,1), x=c(0,1)) + 
          facet_grid(time~model + age) + theme_bw() + labs(x="", y="") +
          scale_fill_gradientn(colors=heat.colors(8)))
dev.off()
