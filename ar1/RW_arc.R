rm(list=ls())
library(ggplot2)
library(data.table)
library(dplyr)
set.seed(123)

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

draws <- 1000
N <- 41
sigma <- 1

ts_df <- data.table(data=c(RW_draws(draws, N, sigma)), time=rep(1:N, draws),
                    draw=rep(1:draws, each=N))


ggplot(ts_df, aes(time, data, group=draw)) + geom_path(alpha=.3)


ts_df <- ts_df %>% group_by(draw) %>% 
                   mutate(lag.10=lag(data, 10), lag.20=lag(data, 20),
                          lag.30=lag(data, 30), lag.40=lag(data, 40)) %>% 
                   as.data.table
ts_df_delta <- subset(ts_df, time==41)


# how do i do this in DT
ts_df_delta[, delta10:= (data - lag.10) / 10]
ts_df_delta[, delta20:= (data - lag.20) / 20]
ts_df_delta[, delta30:= (data - lag.30) / 30]
ts_df_delta[, delta40:= (data - lag.40) / 40]

# and this
dt <- data.table(ARC=c(ts_df_delta$delta10, ts_df_delta$delta20, 
                       ts_df_delta$delta30, ts_df_delta$delta40),
                 draw=rep(1:draws, 4), 
                 time_diff=rep(c("10", "20", "30","40"), each=draws))

ggplot(data = dt, aes(x = as.factor(time_diff), y = ARC)) +
    labs(x = "Time Difference", y = "Annualized Rate of Change")+
    theme_bw() +
    scale_colour_brewer(palette = "Set1")+
    scale_shape_identity() +
    geom_point(aes(shape = 95), size = 7) 

data = c("one"=1, "two"=2, "three"=3)
sapply(data, names)
