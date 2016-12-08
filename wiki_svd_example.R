rm(list=ls())

M <- matrix(0, nrow=4, ncol=5)
M[1,1] <- 1.
M[1,5] <- 2.
M[2,3] <- 3.
M[4,2] <- 2.
M2 <- t(M)

decomp <- svd(M)
decomp$u %*% diag(decomp$d) %*% t(decomp$v)

decomp2 <- svd(M2)
decomp2$u %*% diag(decomp2$d) %*% t(decomp2$v)
