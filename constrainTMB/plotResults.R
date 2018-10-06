rm(list=ls())
library(ggplot2)
library(dplyr)

f_ <- list.files("~/Documents/TMBtests/constrain/results/", full.names=T)

simResultDF <- bind_rows(lapply(f_, read.csv)) %>%
  mutate(constrain=as.logical(constrain))

simResultDF %>%
  ggplot(aes(x=n, y=time.elapsed, color=constrain)) +
  geom_point() +
  coord_trans(x="log", y="log")

simResultDF %>%
  ggplot(aes(x=b0, color=as.factor(constr)))