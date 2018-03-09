library(spatstat)
library(NHPoisson)

data(nztrees)
ppm(nztrees, ~1, Poisson())
# fit the stationary Poisson process to 'nztrees'
# no edge correction needed

lambda <- 1/60 #1 event per minute
time.span <- 60*60*24 #24 hours, with time granularity one second

aux <- simNHP.fun(rep(lambda,time.span))
out<-fitPP.fun(posE=aux$posNH,n=time.span,start=list(b0=0))
