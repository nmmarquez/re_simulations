rm(list=ls())
library(mvtnorm)
library(ggplot2)
library(data.table)
library(dplyr)
library(MASS)
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

time_plot <- function(df, title=""){
    p <- ggplot(df, aes(year_id, y_obs, group=age_group_id, 
                            colour=age_group_id, fill=age_group_id)) +
        geom_path(linetype="dashed", na.rm=T) +
        scale_color_gradientn(colors=rainbow(7)) +
        scale_fill_gradientn(colors=rainbow(7))
    if("y_pred" %in% names(df)){
        p <- p + geom_line(data=df, aes(x=year_id, y=y_pred, 
                                        colour=age_group_id, 
                                        group=age_group_id))
    }
    if("y_lower" %in% names(df)){
        p <- p + geom_ribbon(aes(ymin=y_lower, ymax=y_upper), alpha=0.1, linetype="blank")
    }
    p
}

age_avg <- function(df){
    df_avg <- df[,lapply(.SD, mean), by=year_id]
    df_avg$y_pred <- apply(df_avg[,new_cols, with=F], 1, mean)
    df_avg$y_lower <- apply(df_avg[,new_cols, with=F], 1, quantile, .025)
    df_avg$y_upper <- apply(df_avg[,new_cols, with=F], 1, quantile, .975)
    df_avg
}

age_group_id <- 2:21 # age groups
time_points <- 1990:2015 # number of time points 
forecasting_time_points <- 1990:2040 # forcasted time series
A_ <- length(age_group_id)
T_ <- length(forecasting_time_points)
train_years <- length(time_points)

rho_time <- .98 # correlation over time
rho_age <- .95 # correlation over age
sigma_age <- 1
sigma_time <- 1

Q_age <- Q_ar1(length(age_group_id), sigma_age, rho_age) # precision matrix for age
dimnames(Q_age) <- list(age_group_id, age_group_id)
Q_time <- Q_ar1(length(time_points), sigma_time, rho_time) # precision matrix for time
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
time_plot(df)

# plot average across ages
time_plot(df[,lapply(.SD, mean), by=year_id])

# fit a model
lm1 <- lm(y_obs ~ year_id + age_group_id, data=df)
summary(lm1)

# lets forecast now with a RW on the residuals independent for each age
df$res <- lm1$res
df$y_hat <- lm1$fitted.values
df[, lag.res:=c(NA, res[-.N]), by=age_group_id]
df[, delta.res:=res - lag.res]

rw_df <- df[,list(sigma.rw=sd(delta.res, na.rm=TRUE)), by=age_group_id]

draws <- 1000
df_forecast <- as.data.table(expand.grid(year_id=forecasting_time_points, 
                                         age_group_id=age_group_id))

sim_draws <- do.call(rbind, lapply(rw_df$sigma.rw, function(x)
    RW_draws(draws, T_, x, length(time_points)))) + 
    df_forecast$year_id * lm1$coefficients[2] + lm1$coefficients[1] + 
    df_forecast$age_group_id * lm1$coefficients[3]

new_cols <- paste0("draw", 1:draws)
dimnames(sim_draws) <- list(NULL, new_cols)

df_forecast <- cbind(df_forecast, sim_draws)


# how do I do this in data table?
df_forecast[,y_pred:=rowMeans(.SD), .SDcols=new_cols]
df_forecast$y_lower <- apply(sim_draws, 1, quantile, .025)
df_forecast$y_upper <- apply(sim_draws, 1, quantile, .975)

df_forecast <- as.data.table(left_join(df_forecast, df))

time_plot(df_forecast)

dff_bad <- age_avg(df_forecast)
time_plot(dff_bad)

df_forecast <- as.data.table(expand.grid(year_id=forecasting_time_points, 
                                         age_group_id=age_group_id))
start_points <- 0:(A_ - 1) * T_ + 1
end_points <- 1:A_ * T_
sim_draw_list <- lapply(1:A_, function(i) 
    sim_draws[start_points[i]:end_points[i],])


# how correlated is the age data in the past?
past_data <- matrix(df$y_obs, nrow=train_years)
age_corr <- median(sapply(2:ncol(past_data), function(x) 
    cor(past_data[,x], past_data[,x-1])))


X_age <- t(sapply(sim_draw_list, colSums))
corr_mat_age <- age_corr**abs(outer(0:(A_-1), 0:(A_-1), "-"))
mvdat <- t(mvrnorm(n=draws, mu=0 * 1:A_, Sigma=corr_mat_age, empirical=TRUE))
ranks <- t(apply(mvdat, 1, rank, ties.method="first"))
sorted_X_age <- t(apply(X_age, 1, function(x) sort(x, index.return=TRUE)$ix))
sorted_age_sims <- lapply(1:A_, function(x) 
    sim_draw_list[[x]][,sorted_X_age[x,]])

age_time_corr_X <- lapply(1:A_, function(x)
    sorted_age_sims[[x]][,ranks[x,]])

sim_draws_corr <- do.call(rbind, age_time_corr_X)
dimnames(sim_draws_corr) <- list(NULL, new_cols)
df_forecast <- cbind(df_forecast, sim_draws_corr)


# how do I do this in data table?
df_forecast[,y_pred:=rowMeans(.SD), .SDcols=new_cols]
df_forecast$y_lower <- apply(sim_draws, 1, quantile, .025)
df_forecast$y_upper <- apply(sim_draws, 1, quantile, .975)

df_forecast <- as.data.table(left_join(df_forecast, df))

# show past data here too
time_plot(df_forecast)

dff <- age_avg(df_forecast)

time_plot(dff)

