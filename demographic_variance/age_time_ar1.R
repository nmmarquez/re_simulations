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
    RW_draws(draws, length(forecasting_time_points), x, length(time_points)))) + 
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

# show past data here too
time_plot(df_forecast)
    
dff <- df_forecast[,lapply(.SD, mean), by=year_id]
dff$y_pred <- apply(dff[,new_cols, with=F], 1, mean)
dff$y_lower <- apply(dff[,new_cols, with=F], 1, quantile, .025)
dff$y_upper <- apply(dff[,new_cols, with=F], 1, quantile, .975)

time_plot(dff)

time_plot(subset(df_forecast, age_group_id == 11))


Q_af <- Q_ar1(length(age_group_id), sigma_Q, rho_age) 
dimnames(Q_age) <- list(age_group_id, age_group_id)
Q_tf <- Q_ar1(length(forecasting_time_points), sigma_Q, rho_time) 
dimnames(Q_time) <- list(time_points, time_points)

Qf <- kronecker(Q_af, Q_tf, make.dimnames=TRUE) # joint precision matrix


# observe values which have correlated error
true_process <- lm1$coefficients[1] + lm1$coefficients[2] * df_forecast$year_id + 
    lm1$coefficients[3] * df_forecast$age_group_id + 
    t(rmvnorm(draws, sigma=solve(Qf)))

df_true <- as.data.table(expand.grid(year_id=forecasting_time_points, 
                                     age_group_id=age_group_id))

for(k in 1:length(new_cols)){
    df_true[,(new_cols[k]):=true_process[,k],]
}

df_true$ref <- lm1$coefficients[1] + lm1$coefficients[2] * df_forecast$year_id + 
  lm1$coefficients[3] * df_forecast$age_group_id

df_adj <- subset(df_true, year_id == 2015)
for(k in 1:length(new_cols)){
    df_adj[,(new_cols[k]):=df_adj[[new_cols[k]]] - df_adj[["ref"]],]
}

df_adj$year_id <- NULL
df_adj$ref <- NULL
df_true_merged <- as.data.table(left_join(df_true, df_adj, by="age_group_id"))

for(k in 1:length(new_cols)){
    df_true_merged[,(new_cols[k]):=df_true_merged[[paste0(new_cols[k], ".x")]] - df_true_merged[[paste0(new_cols[k], ".y")]],]
}

# how do I do this in data table?
df_true_merged$y_pred <- apply(df_true_merged[,new_cols,with=F], 1, mean)
df_true_merged$y_lower <- apply(df_true_merged[,new_cols,with=F], 1, quantile, .025)
df_true_merged$y_upper <- apply(df_true_merged[,new_cols,with=F], 1, quantile, .975)
df_true_merged[year_id < 2015, y_lower:=y_pred,]
df_true_merged[year_id < 2015, y_upper:=y_pred,]

df_true_merged <- as.data.table(left_join(df_true_merged, df))
time_plot(df_true_merged)

time_plot(subset(df_true_merged, age_group_id == 10))

dfft <- df_true_merged[,lapply(.SD, mean), by=year_id]
dfft$y_pred <- apply(dfft[,new_cols, with=F], 1, mean)
dfft$y_lower <- apply(dfft[,new_cols, with=F], 1, quantile, .025)
dfft$y_upper <- apply(dfft[,new_cols, with=F], 1, quantile, .975)
dfft[year_id < 2015, y_lower:=y_pred,]
dfft[year_id < 2015, y_upper:=y_pred,]

time_plot(dfft)
time_plot(dff)
