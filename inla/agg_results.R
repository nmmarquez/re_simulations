rm(list=ls())
pacman::p_load(ggplot2, data.table)

save_dir <- "/home/j/temp/nmarquez/sim2dresults/"
metaf <- paste0(save_dir, grep("_meta.Rda", list.files(save_dir), value=TRUE))

DT <- data.table()

for(m in metaf){
    load(m)
    names(MetaList) <- c("rho", "kappa", "tau", "sigma", "range", "N",
                         "inrmsediff", "outrmsediff", "inmaddiff", "outmaddiff", 
                         "inla.time", "tmb.time")
    MetaList[["inla.time"]] <- as.numeric(MetaList$inla.time, units="mins")
    MetaList[["tmb.time"]] <- as.numeric(MetaList$tmb.time, units="mins")
    DT <- rbindlist(list(DT, as.data.table(MetaList)), fill=TRUE)
}

DT[,diff.time:=inla.time-tmb.time]
DT[,ratio.time:= as.numeric(inla.time) / as.numeric(tmb.time)]
DT[,N:=as.factor(N)]

jpeg("~/Desktop/inla_tmb_time_diff.jpg", width = 960, height = 960)
ggplot(DT, aes(x=diff.time, fill=N, group=N)) + geom_density(alpha=.4) +
    geom_vline(xintercept=0) + xlab("INLA time - TMB time(minutes)")
dev.off()

jpeg("~/Desktop/inla_tmb_time_ratio.jpg", width = 960, height = 960)
ggplot(DT, aes(x=ratio.time, fill=N, group=N)) +
    geom_density(alpha=.4) + geom_vline(xintercept=1) +
    xlab("INLA time / TMB time(minutes)") 
dev.off()

jpeg("~/Desktop/inla_tmb_time_ratio_kappa.jpg", width = 960, height = 960)
ggplot(DT[kappa < 5 | kappa > 6], aes(x=ratio.time, fill=N, group=N)) +
    geom_density(alpha=.4) + geom_vline(xintercept=1) +
    xlab("INLA time / TMB time(minutes)") + facet_wrap(~kappa)
dev.off()

jpeg("~/Desktop/rmse_diff.jpg", width = 960, height = 960)
ggplot(DT[outrmsediff > -.1], aes(x=outrmsediff, fill=N, group=N)) +
    geom_density(alpha=.4) + xlab("") +
    labs(title="RMSE Difference(INLA-TMB)")
dev.off()

jpeg("~/Desktop/maddiff.jpg", width = 960, height = 960)
ggplot(DT[outmaddiff > -.1], aes(x=outmaddiff, fill=N, group=N)) +
    geom_density(alpha=.4) + xlab("") +
    labs(title="MAD Difference(INLA-TMB)")
dev.off()

jpeg("~/Desktop/time_diff_point.jpg", width = 960, height = 960)
ggplot(DT, aes(x=inla.time, y=tmb.time)) + geom_point(alpha=.4) + 
    facet_wrap(~N) + geom_abline(intercept=0, slope=1) +
    labs(title="Time Comparison By Sample Size")
dev.off()
