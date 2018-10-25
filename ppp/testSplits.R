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
set.seed(1234)

US.df$COUNTRY <- 1

boundSPDF <- gUnaryUnion(US.df, id=US.df$COUNTRY) %>%
    SpatialPolygonsDataFrame(
        data.frame(
            COUNTRY=row.names(.),
            row.names=row.names(.)))

voronoi.custom <- function (x, bounding.polygon, ...){
    if (!inherits(x, "SpatialPoints")) {
        stop("Must pass a SpatialPoints* object to voronoi.polygons.")
    }
    crds <- coordinates(x)
    z <- deldir::deldir(crds[, 1], crds[, 2], ...)
    w <- deldir::tile.list(z)
    polys <- vector(mode = "list", length = length(w))
    for (i in seq(along = polys)) {
        pcrds <- cbind(w[[i]]$x, w[[i]]$y)
        pcrds <- rbind(pcrds, pcrds[1, ])
        polys[[i]] <- Polygons(list(Polygon(pcrds)), ID = as.character(i))
    }
    SP <- SpatialPolygons(polys, proj4string = CRS(proj4string(x)))
    voronoi <- SpatialPolygonsDataFrame(
        SP, 
        data = data.frame(
            x = crds[,1], 
            y = crds[, 2], 
            area = sapply(slot(SP, "polygons"),                                                                                            slot, "area"), row.names = sapply(slot(SP, "polygons"), 
                                                                                                                                         slot, "ID")))
    if (!missing(bounding.polygon)) {
        bounding.polygon <- gUnion(bounding.polygon, bounding.polygon)
        voronoi.clipped <- gIntersection(voronoi, bounding.polygon, 
                                         byid = TRUE, id = row.names(voronoi))
        df <- data.frame(voronoi)
        df$area <- sapply(slot(voronoi.clipped, "polygons"), slot, "area")
        voronoi <- SpatialPolygonsDataFrame(voronoi.clipped, df)
    }
    voronoi
}

# Create n sample
N <- c(10, 50, 100, 200)
names(N) <- N

if(!file.exists("~/Documents/re_simulations/ppp/polysList.Rds")){
    polysList <- lapply(N, function(i){
        tests <- sapply(1:1000, function(j){
            sampleSPDF <- spsample(boundSPDF, i, "random", iter=10000)
            voronoi.custom(
                sampleSPDF, 
                boundSPDF, 
                rw=c(t(boundSPDF@bbox)),
                eps=1e-06)
        })
        testMaxRatios <- sapply(tests, function(x) {
            max(outer(x$area, x$area, `/`))})
        tests[[which(testMaxRatios == min(testMaxRatios))]]
    })
    
    saveRDS(polysList, file="~/Documents/re_simulations/ppp/polysList.Rds")

}

polysList <- read_rds("~/Documents/re_simulations/ppp/polysList.Rds")