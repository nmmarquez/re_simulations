library(abind)

N <- 7 # number of observations
M <- 2 # number of covariates
D <- 10 # number of draws

# simulate the covariates in 3-D array
covs <- array(rnorm(N * M * D), dim =c(N, M, D)) 
covs[,1,] <- 1 # make sure the first cov is always 1 (the intercept)
dim(covs)

# simulate some betas
betas <- matrix(c(rnorm(D,3), rnorm(D, -1)), nrow=D)
# make teh arrays conformable
betas_expand <- array(rep(c(t(betas)), each=N), dim=c(N, M, D))

# they conform!
dim(covs)
dim(betas_expand)

# sum over the 2nd axis
apply(covs * betas_expand, c(1,3), sum)

t(betas)
betas_expand2 <- array(t(betas), dim=(c(1, M, D)))
apply(covs * betas_expand2, c(1,3), sum)
