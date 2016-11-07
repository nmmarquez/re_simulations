rm(list=ls())
set.seed(124)
library(gee)
library(lme4)
library(mvtnorm)
library(ggplot2)

options(lmerControl=list(check.nobs.vs.rankZ = "warning",
                         check.nobs.vs.nlev = "warning",
                         check.nobs.vs.nRE = "warning",
                         check.nlev.gtreq.5 = "warning",
                         check.nlev.gtr.1 = "warning"))

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

# summaries of each model type
get_betas <- function(model){
    if(isS4(model)){
        b <- model@beta
    }
    else{
        b <- model$coefficients
    }
    b
}

run_model <- function(fun=c(lm, lmer, gee), df=data,
                      ff=c("y ~ x1", "y ~ x1 + (1|location)", "y ~ x1"),
                      model_names=c("lm", "lmer", "gee"),
                      iterations=100, res="unstructured"){
    sample_sizes <- as.integer(seq(0, N , N/iterations))
    # make sure samples are greater than 10
    location <- df$location
    sample_sizes <- sample_sizes[sample_sizes >= 10]
    samples <- lapply(sample_sizes, function(x) sample(1:N, size=x))
    M <- length(fun)
    models <- lapply(1:M, function(x) t(sapply(samples, function(y)
        get_betas(fun[[x]](ff[[x]], id=location, data=df[y,], corstr=res)))))
    betas <- paste0("beta", 0:(ncol(models[[1]]) - 1))
    results <- data.frame(do.call(rbind, models))
    names(results) <- betas
    results$sample_size <- rep(sample_sizes, M)
    results$model <- rep(model_names, each=length(sample_sizes))
    results
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
data <- data.frame(y, x1, location=site_number_sample)
data <- data[order(data$location),]


beta_results <- run_model()

# blackline is the true value
ggplot(data = beta_results, aes(x=sample_size, y=beta0, colour=model)) +       
    geom_line() + geom_abline(aes(slope=0, intercept=b0)) +
    ggtitle("B0 estimates for Correct Residual Var Structure of LMER")
ggplot(data = beta_results, aes(x=sample_size, y=beta1, colour=model)) +       
    geom_line() + geom_abline(aes(slope=0, intercept=b1)) + 
    ggtitle("B1 estimates for Correct Residual Var Structure of LMER")

# now use the wrong covariance structure in the random effects model
ar1_rho <- .99
ar1_sigma <- .1
ar1_re_eff <- sim_ar1(site_number, ar1_rho, ar1_sigma)
qplot(1:site_number, ar1_re_eff, main="AR1 Correlation Over Time")
ar1_site_eff <- ar1_re_eff[site_number_sample]
b2 <- 1.5

# generate the response variables
y2 <- b0 + b1 * x1 + b2 * site_number_sample + ar1_site_eff + rnorm(N, sd=.2)
data2 <- data.frame(y, y2, x1, location=site_number_sample)
data2 <- data2[order(data2$location),]

beta_results2 <- run_model(df=data2,
                          ff=c("y2 ~ x1 + location", 
                               "y2 ~ x1 + location + (1|location)", 
                               "y2 ~ x1 + location"))

# summaries of each model type
ggplot(data = beta_results2, aes(x=sample_size, y=beta0, colour=model)) +       
    geom_line() + geom_abline(aes(slope=0, intercept=b0)) +
    ggtitle("B0 estimates for Incorrect Residual Var Structure of LMER")
ggplot(data = beta_results2, aes(x=sample_size, y=beta1, colour=model)) +       
    geom_line() + geom_abline(aes(slope=0, intercept=b1)) + 
    ggtitle("B1 estimates for Incorrect Residual Var Structure of LMER")
ggplot(data = beta_results2, aes(x=sample_size, y=beta2, colour=model)) +       
    geom_line() + geom_abline(aes(slope=0, intercept=b2)) + 
    ggtitle("B2 estimates for Incorrect Residual Var Structure of LMER")
