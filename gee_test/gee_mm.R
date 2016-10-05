rm(list=ls())
set.seed(123)
library(gee)
library(lme4)
library(mvtnorm)

# function for generating an AR1 precision and or var-covar matrix
ar1_matrix <- function(N, rho, sigma, vcov=FALSE){
    mat <- matrix(data=0, nrow=N, ncol=N)
    for(i in 2:N){
        mat[i,i] <- (1. + rho**2) / sigma**2
        mat[i-1,i] <- -rho / sigma**2
        mat[i,i-1] <- -rho / sigma**2
    }
    for(i in c(1,N)){
        mat[i,i] <- 1. / sigma**2
    }
    if(vcov){
        mat <- solve(mat)
    }
    mat
}

# simulate an ar1 process
sim_ar1 <- function(N, rho, sigma){
    vcov_mat <- ar1_matrix(N, rho, sigma, vcov=TRUE)
    as.vector(rmvnorm(1, mean=rep(0, N), sigma=vcov_mat))
}

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

run_model <- function(iterations=200, res="unstructured"){
    sample_sizes <- seq(0, N , N/iterations)
    var <- 
    # make sure samples are greater than 20
    sample_sizes <- sample_sizes[sample_sizes >= 20]
    samples <- lapply(sample_sizes, function(x) sample(1:N, size=x))
    m1 <- t(sapply(samples, function(x) lm(y ~ x1, data=df[x,])$coefficients))
    m2 <- t(sapply(samples, function(x)
        lmer(y ~ x1 + (1|location), data=df[x,])@beta))
    m3 <- t(sapply(samples, function(x)
        gee(y ~ x1, id=location, data=df[x,], corstr=res)$coefficients))
    return(None)
}
pop_samp <- sample(1:N, site_number)
m1 <- lm(y ~ x1)
m2 <- lmer(y ~ x1 + (1|location), data=df)
m3 <- gee(y ~ x1, id=location, data=df, corstr="unstructured")

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
