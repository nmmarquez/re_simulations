rm(list=ls())
library(tidyverse)
library(sf)

setwd("~/Documents/re_simulations/ppp/data/")

# everything in here is a cbh
allDF <- read_rds("./prepared_cbh_2018_09_04_11_49_25.RDS")
head(allDF)

# looking at this group dominican republic("DOM") looks like a good candidate
allDF %>%
    mutate(has_coords=is.na(longnum)) %>%
    group_by(country, has_coords) %>%
    summarize(n=n()) %>%
    filter(n() == 2)

DF <- filter(allDF, country == "DOM") %>%
    as_tibble %>%
    mutate(GAUL_CODE=as.numeric(location_code))

cat(paste0("Number of Observations: ", sum(DF$N), "\n"))
cat(paste0("Number of Deaths: ", round(sum(DF$died)), "\n"))
cat(paste0("Proportion: ", round(sum(DF$died)/sum(DF$N), 4), "\n"))

admin1Codes <- DF %>%
    filter(shapefile == "admin2013_1") %>%
    select(location_code) %>%
    unique %>%
    unlist %>%
    unname

# row counts
table(DF$shapefile)

# weighted counts
DF %>%
    group_by(shapefile) %>%
    summarize(personCount=sum(N))

DF %>%
    select(nid, shapefile) %>%
    rename(Nid=nid) %>%
    unique %>%
    left_join(read_csv("./id.csv"), by="Nid") %>%
    select(Title, shapefile)

# temporal breakdown: I definitely dont get this
DF %>%
    group_by(shapefile, year, svyyr, source) %>%
    summarize(personCount=sum(N))

shapeADMIN <- st_read("admin2013_1.shp")
shapeADMINDR <- shapeADMIN[shapeADMIN$GAUL_CODE %in% admin1Codes,]

shapeFiles <- grep("DOM", unique(DF$shapefile), value=T)
shapeList <- c(list(shapeADMINDR), lapply(paste0(shapeFiles, ".shp"), st_read))
names(shapeList) <- c("admin2013_1", shapeFiles)

dataCountList <- lapply(names(shapeList), function(spn){
    countDF <- DF %>%
        group_by(shapefile, GAUL_CODE) %>%
        summarize(obsCount=sum(N)) %>%
        filter(shapefile == spn) %>%
        ungroup %>%
        select(-shapefile)
    left_join(shapeList[[spn]], countDF, by="GAUL_CODE") %>%
        mutate(set=spn) %>%
        select(set, GAUL_CODE, obsCount, geometry)
})

names(dataCountList) <- names(shapeList)

# for right now lets just plot the data that we do have
pointSFDF <- DF %>% 
    filter(!is.na(longnum)) %>%
    st_as_sf(
        coords = c("longnum", "latnum"), 
        crs = "+proj=longlat +datum=WGS84 +no_defs") %>%
    st_combine

# Beautiful!
ggplot(DF, aes(x=longnum, y=latnum)) +
    geom_point() + 
    theme_classic()

# No :(
ggplot(pointSFDF) +
    geom_sf(data=shapeADMINDR) +
    geom_sf(size=.2, alpha=.3) +
    theme_classic()

do.call(rbind, dataCountList) %>%
    ggplot(aes(fill=obsCount)) +
    geom_sf() +
    scale_fill_distiller(palette="Spectral") +
    theme_classic() +
    facet_wrap(~set, nrow=2)

do.call(rbind, dataCountList[1:3]) %>%
    ggplot(aes(fill=obsCount)) +
    geom_sf() +
    scale_fill_distiller(palette="Spectral") +
    theme_classic() +
    facet_wrap(~set, nrow=2)
