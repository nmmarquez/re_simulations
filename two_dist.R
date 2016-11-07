# sample code for examing bimodal distributions
rm(list=ls())
set.seed(123)

create_distribution <- function(N, p, mean1=0, mean2=0, sd1=1, sd2=1){
    funcs <- list((function (N) rnorm(N, mean1, sd1)), 
                   (function (N) rnorm(N, mean2, sd2)))
    draws <- table(rbinom(N, 1, p) + 1)
    c(funcs[[1]](draws[1]), funcs[[2]](draws[2]))
}

N <- 10**6
samples <- create_distribution(N, .5, mean1=0, sd1=2, mean2=3, sd2=2)
mean(samples)
sd(samples)
hist(samples)
