rm(list = ls())
library(httr)
library(jsonlite)
library(dplyr)

orskey <- read.table("~/.openrouteservicetoken", stringsAsFactors=F)$V1
ghkey <- read.table("~/.graphhopperkey", stringsAsFactors=F)$V1

# test results against open street map interface
# https://www.openstreetmap.org/#map=14/47.6385/-122.3331

# The geocoding api through nominatim docs no key necessary
# https://wiki.openstreetmap.org/wiki/Nominatim
paste0("https://nominatim.openstreetmap.org/",
       "search?city=los+angeles&state=california&format=json") %>%
    GET %>% prettify %>% fromJSON

# docs for open route services which you need an api key for
# https://go.openrouteservice.org/documentation/
paste0("https://api.openrouteservice.org/directions?",
       "api_key=", orskey, "&profile=cycling-regular&units=km",
       "&coordinates=-122.3271%2C47.6226%7C-122.3080%2C47.6543") %>%
    GET %>% prettify %>% fromJSON %>% str

# docs for graph hopper which you need api key for
# https://graphhopper.com/api/1/docs/routing/
# time is in milliseconds by default why???
paste0("https://graphhopper.com/api/1/route?",
       "point=47.6226%2C-122.3271&point=47.6543%2C-122.3080",
       "&vehicle=bike&locale=en&key=", ghkey) %>%
    GET %>% prettify %>% fromJSON %>% str
