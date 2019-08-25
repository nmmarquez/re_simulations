rm(list=ls())
library(ar.matrix)
library(TMB)
library(sp)
library(leaflet)
library(Matrix)

set.seed(123)
US.df@data$zeta <- c(r.lCAR(1, graph=US.graph, sigma=10, rho=.9999))
US.df@data$iid <- rnorm(nrow(US.df@data), 0, .3)
US.df@data$data <- US.df@data$iid + US.df@data$zeta 

pal <- colorNumeric(palette="YlGnBu", domain=US.df@data$data)

# see map
map1<-leaflet() %>%
    addProviderTiles("CartoDB.Positron") %>%
    addPolygons(data=US.df, fillColor=~pal(data), color="#b2aeae",
                fillOpacity=0.7, weight=0.3, smoothFactor=0.2) %>%
    addLegend("bottomright", pal=pal, values=US.df$data, title="", opacity=1)

map1

modelRun <- function(y, verbose=T, rho=1){
    model <- "icar"
    compile(paste0(model, ".cpp"))
    dyn.load(dynlib(model))
    Data <- list(
        yobs=y, Wstar=Matrix(diag(rowSums(US.graph)) - US.graph, sparse=T))
    dim(Data$yobs) <- length(y)
    Params <- list(log_sigma=0, log_sigmasp=0, zeta=rep(0, length(y)), rho=rho)
    
    Map <- list(rho=factor(NA))
    
    Obj <- MakeADFun(data=Data, parameters=Params, DLL=model, random="zeta",
                     silent=!verbose, map=Map)

    Opt <- nlminb(start=Obj$par, objective=Obj$fn, 
                  gradient=Obj$gr,
                  control=list(eval.max=1e6, iter.max=1e6))

    sdrep <- sdreport(Obj, getJointPrecision=TRUE)
    dyn.unload(dynlib(model))

    list(obj=Obj, opt=Opt, sd=sdrep)
}


test <- modelRun(US.df$data, rho=1.)
