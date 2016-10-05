rm(list=ls())
library(ggplot2)

ar1_sim <- function(N, rho, sigma, trend=0){
    obs <- 1:N * trend
    for(i in 2:N){
        res_past <- obs[i-1] - ((i-1) * trend)
        obs[i] <- obs[i] + rnorm(1, mean=rho * res_past, sd=sigma)
    }
    obs
}

time_point_var <- function(N_obs, N_forecast, rho, sigma){
    vars <- rep(0, N_obs + N_forecast)
    for(i in (N_obs + 1):(N_obs+N_forecast)){
        vars[i] <- vars[i-1] + sigma**2 * rho**(2 * (i-N_obs-1))
    }
    vars
}

ar1_forecast <- function(N_obs, N_forecast, rho, sigma, trend=0){
    t <- 1:(N_obs + N_forecast)
    obs <- ar1_sim(N_obs, rho, sigma, trend)
    obs_trend <- ifelse(trend !=0, (obs[N_obs] - obs[1]) / N_obs, obs)
    y <- c(obs, rep(obs[N_obs], N_forecast) + 1:N_forecast * obs_trend)
    sds <- time_point_var(N_obs, N_forecast, rho, sigma)**.5
    data.frame(t=t, y=y, lower_bound=y - 1.96 * sds, upper_bound=y + 1.96 * sds)
}

N_obs <- 20
N_forecast <- 20
sigma <- 1
rho <- .99
trend <- 0
df <- ar1_forecast(N_obs, N_forecast, rho, sigma, trend)
ggplot(df, aes(x=t, y=y)) + geom_line() + 
    geom_ribbon(aes(ymin=lower_bound, ymax=upper_bound), alpha=.4)
