library(mvtnorm)

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

sim_ar1 <- function(N, rho, sigma){
    vcov_mat <- ar1_matrix(N, rho, sigma, vcov=TRUE)
    as.vector(rmvnorm(1, mean=rep(0, N), sigma=vcov_mat))
}
