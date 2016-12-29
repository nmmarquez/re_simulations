rm(list=ls())
library(mvtnorm)
library(ggplot2)
library(data.table)
library(dplyr)
set.seed(123)

Q_ar1 <- function(N, sigma, rho){
    Q <- matrix(0, nrow=N, ncol=N)
    Q[1,1] <- 1.
    for(i in 2:N){
        Q[i,i] <- 1 + rho**2
        Q[i-1,i] <- -1 * rho
        Q[i,i-1] <- -1 * rho
    }
    Q[N,N] <- 1.
    (1 / sigma**2) * Q
}

RW <- function(N, sigma, burn_in=1){
    time_series <- rnorm(N, mean=0, sd=sigma)
    for(i in (burn_in + 1):N){
        time_series[i] <- rnorm(1, mean=time_series[i-1], sd=sigma)
    }
    return(time_series)
}

RW_draws <- function(draws, N, sigma, burn_in=1){
    sapply(1:draws, function(x) RW(N, sigma, burn_in))
}

age_group_id <- 2:21 # age groups
time_points <- 1990:2015 # number of time points 
forecasting_time_points <- 1990:2040 # forcasted time series
train_years <- max(time_points)

rho_time <- .98 # correlation over age
rho_age <- .95 # correlation over
sigma_Q <- 1


Q_age <- Q_ar1(length(age_group_id), sigma_Q, rho_age) # precision matrix for age
dimnames(Q_age) <- list(age_group_id, age_group_id)
Q_time <- Q_ar1(length(time_points), sigma_Q, rho_time) # precision matrix for time
dimnames(Q_time) <- list(time_points, time_points)

Q <- kronecker(Q_age, Q_time, make.dimnames=TRUE) # joint precision matrix
df <- as.data.table(expand.grid(year_id=time_points, age_group_id=age_group_id))

B_time <- 2 # beta on time
B_age <- 2 # beta on age
B0 <- -1990 * 2 # intercept

# observe values which have correlated error
df$y_obs <- B0 + B_time * df$year_id + B_age * df$age_group_id + 
    c(rmvnorm(1, sigma=solve(Q)))

# plot age specific
ggplot(df, aes(year_id, y_obs, group=age_group_id, colour=age_group_id)) +
    geom_path(alpha = 0.5) + scale_color_gradientn(colors=rainbow(7))

# plot average across ages
ggplot(df[,mean(y_obs), by=year_id], aes(year_id, V1)) + geom_path(alpha=0.5)

# fit a model
lm1 <- lm(y_obs ~ year_id + age_group_id, data=df)
summary(lm1)

# lets forecast now with a RW on the residuals independent for each age
df$res <- lm1$res
df$y_hat <- lm1$fitted.values
df[, time.lag.res:=c(NA, res[-.N]), by=age_group_id]
df[, time.delta.res:=res - time.lag.res]
df <- df[order(year_id, age_group_id)]
df[, age.lag.res:=c(NA, res[-.N]), by=year_id]
df[, age.delta.res:=res - age.lag.res]

df <- df[order(age_group_id, year_id)]

rwt_df <- df[,list(sigma.rw=sd(time.delta.res, na.rm=TRUE)), by=age_group_id]
rwa_df <- df[,list(sigma.rw=sd(age.delta.res, na.rm=TRUE))]
rwa_df
