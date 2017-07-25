rm(list=ls())
library(ggplot2)

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
    DT <- rbind(DT, as.data.table(MetaList))
}

DT[,diff.time:=inla.time-tmb.time]
DT[,ratio.time:= as.numeric(inla.time) / as.numeric(tmb.time)]
DT[,N:=as.factor(N)]

ggplot(DT, aes(x=diff.time, fill=N, group=N)) + geom_density(alpha=.4) + 
    geom_vline(xintercept=0) + xlab("INLA time - TMB time(minutes)")
ggplot(DT, aes(x=ratio.time, fill=N, group=N)) + geom_density(alpha=.4) + 
    geom_vline(xintercept=1) + xlab("INLA time / TMB time(minutes)") +
    facet_wrap(~sigma)

ggplot(DT[outrmsediff > -.1], aes(x=outrmsediff, fill=N, group=N)) + 
    geom_density(alpha=.4) + xlab("")
ggplot(DT[outmaddiff > -.1], aes(x=outrmsediff, fill=N, group=N)) + 
    geom_density(alpha=.4) + xlab("")
