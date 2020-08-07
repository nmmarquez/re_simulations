rm(list=ls())
library(tidyverse)
library(tigris)
library(sf)

state_vec <- c(state.abb, "DC")

# download all the states county data
county_df <- do.call(
    sf:::rbind.sf, lapply(state_vec, counties, class = "sf")) %>%
    select(COUNTYID = GEOID, COUNTY = NAME, geometry)

# download all the states place data
place_df <- do.call(
    sf:::rbind.sf, lapply(state_vec, places, class = "sf")) %>%
    select(PLACEID = GEOID, PLACE = NAME, STATEFP, geometry)

# overlap the two data sets and calculate the area of the place that falls 
# within a particular county over the entire area of a space
# its important to add some buffer so we ignore boundaries that touch
ov_df <- st_intersection(county_df, st_buffer(place_df, 1e-5)) %>%
    # calculate the area of each place county overlap
    mutate(Area = as.numeric(st_area(.))) %>%
    # find out how much each place falls into each county
    group_by(PLACE, STATEFP) %>%
    mutate(pcounty = Area/sum(Area)) %>%
    ungroup()

# here are the locations that have a "significant" portion of their area in 
# multiple counties
ov_df %>%
    filter(pcounty > .05 & pcounty < .95) %>%
    arrange(PLACE, STATEFP) %>%
    as_tibble() %>%
    select(PLACE, STATEFP, COUNTY, pcounty)

# lets save a csv of a places majority county
ov_df %>%
    as_tibble() %>%
    group_by(PLACE, STATEFP) %>%
    filter(pcounty == max(pcounty)) %>%
    select(PLACE, COUNTY, PLACEID, COUNTYID, pcounty, STATEFP) %>%
    ungroup() %>%
    write_csv("majority_county_place.csv")

# and this one will save all counties associated with a place
ov_df %>%
    as_tibble() %>%
    select(PLACE, COUNTY, PLACEID, COUNTYID, pcounty, STATEFP) %>%
    write_csv("full_county_place.csv")
