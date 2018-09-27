rm(list=ls())
pacman::p_load(INLA, ggplot2, data.table, lattice, arm, dplyr, TMB, ar.matrix)
set.seed(124)

# what is the shape that we are dealing with?
plot(US.df)

randomSPDF <- spsample(US.df, 800, "random")

points(randomSPDF, pch=20, col="red", cex=.5)
