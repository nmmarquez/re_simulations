rm(list=ls())
set.seed(123)
source("./ar1.R")
library(gee)
library(lme4)

# sim data with known residual structure 
N <- 1000
site_number <- 100
re_sigma <- 2
re_eff <- rnorm(site_number, sd=re_sigma)
site_number_sample <- sample(1:site_number, size=N, replace=TRUE)
site_eff <- re_eff[site_number_sample]
x1 <- rnorm(N, mean=6, sd=.4)

# beta coefficients
b0 <- 2
b1 <- -3

# generate the response variables
y <- b0 + b1 * x1 + site_eff + rnorm(N, sd=.2)

# put data into data frame for gee
df <- data.frame(y, x1, location=site_number_sample)
df <- df[order(df$location),]

# summaries of each model type
pop_samp <- sample(1:N, site_number)
summary(lm(y ~ x1, data=df[pop_samp,]))
summary(lmer(y ~ x1 + (1|location), data=df[pop_samp,]))
summary(gee(y ~ x1, id=location, data=df[pop_samp,]))

# now use the wrong covariance structure in the random effects model
ar1_rho <- .99
ar1_sigma <- .1
ar1_re_eff <- sim_ar1(site_number, ar1_rho, ar1_sigma)
plot(1:site_number, ar1_re_eff)
ar1_site_eff <- ar1_re_eff[site_number_sample]
b3 <- 1.5

# generate the response variables
y2 <- b0 + b1 * x1 + b3 * site_number_sample + ar1_site_eff + rnorm(N, sd=.2)
df2 <- data.frame(y, y2, x1, location=site_number_sample)
df2 <- df2[order(df2$location),]

# summaries of each model type
summary(lm(y2 ~ x1 + location, data=df2[pop_samp,]))
summary(lmer(y2 ~ x1 + location + (1|location), data=df2[pop_samp,]))
summary(gee(y2 ~ x1 + location, id=location, data=df2, subset=pop_samp))
