rm(list=ls())
library(MASS)
library(tidyverse)
set.seed(123)

M <- 1000
b <- log10(rgamma(M, 6, 3))
a <- log10((rweibull(M, 8, 2)))
clrs <- colorRampPalette(c('white','blue','yellow','red','darkred'))
N <- 250

## ggplot version 
bind_rows(lapply(seq(.09, .04, by=-.01), function(i){
    dsub <- kde2d(a, b, n=N, h=i)
    tibble(z=c(dsub$z), x=rep(dsub$x, N), y=rep(dsub$y, each=N)) %>%
        mutate(bandwidth=i)})) %>%
    mutate(bandwidth=factor(bandwidth, levels=seq(.09, .04, by=-.01))) %>%
    group_by(bandwidth) %>%
    mutate(z = (z-mean(z))/sd(z)) %>%
    ggplot(aes(x=x, y=y, fill=z)) +
    geom_raster() +
    scale_fill_gradientn(colors=clrs(27)) +
    theme_classic() +
    facet_wrap(~bandwidth) +
    guides(fill=FALSE) +
    ggtitle("KDE with varying Bandwidth")
