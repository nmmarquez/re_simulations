rm(list=ls())
library(tidyverse)
library(sp)
library(rgeos)
library(deldir)
library(measurements)
library(mapproj)
library(spatstat)
library(maptools) 
library(SDraw)
library(ar.matrix)


US.df$COUNTRY <- 1

boundSPDF <- gUnaryUnion(US.df, id=US.df$COUNTRY) %>%
    SpatialPolygonsDataFrame(
        data.frame(
            COUNTRY=row.names(.),
            row.names=row.names(.)))

# Create n sample
sampleSPDF <- spsample(boundSPDF, 3, "random")
samplePoly <- voronoi.polygons(sampleSPDF)#, boundSPDF)
plot(samplePoly)
#testSPDF <- gIntersection(boundSPDF, samplePoly)
#plot(testSPDF)