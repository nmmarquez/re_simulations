set.seed(123)
rm(list=ls())
library(forecast)

# simulate ar1 error function
ar1_sim <- function(N, rho, sigma, trend=0){
    obs <- 1:N * trend
    err <- rnorm(N, 0, sigma)
    for(i in 2:N){
        err[i] <- rnorm(1, err[i-1] * rho, sigma) 
    }
    obs + err
}

N <- 1000 # Number of observations
X <- 1:N # Observation Events
rho <- .95 # Corrlation over time
sigma_ar1 <- 4. # correlated variation 
sigma_epsilon <- .1 # random error
trend <- 1. # secular trend
# simulate data and turn into time series object
Y <- ts(ar1_sim(N, rho, sigma_ar1, trend) + rnorm(N, sd=sigma_epsilon))
plot(X, Y, type="l")

# process is non stationy and probably wont fit due to how forecast calcs ar1
fit <- Arima(Y, order=c(1,0,0), include.drift=FALSE, include.constant=FALSE)

# process has an AR of close to 1 may fit but with not great forecats
fit <- Arima(Y, order=c(1,0,0), include.drift=FALSE, include.constant=TRUE)
summary(fit)
plot(forecast(fit, h=200))

# rw process with drift
fit <- Arima(Y, order=c(0,1,0), include.drift=TRUE, include.constant=TRUE)
summary(fit)
plot(forecast(fit, h=200))

# capturing the true process which after accounting for trend is stationary
fit <- Arima(Y, order=c(1,0,0), include.drift=TRUE, include.constant=TRUE)
summary(fit)
plot(forecast(fit, h=200))

