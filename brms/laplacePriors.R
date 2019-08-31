rm(list=ls())

library(brms)
library(dplyr)
library(LaplacesDemon)

set.seed(123)
N <- 100
DF <- data.frame(x = sample(c(FALSE, TRUE), N, replace = TRUE)) %>%
    mutate(w=rnorm(N, c(.3, .7)[x+1], .1))
test <- prior(horseshoe(df=1), class = "b")
m1 <- brm(w ~ 1 + x, data=DF, prior = test)

DF2 <- data.frame(x = sample(c(FALSE, TRUE), N, replace = TRUE)) %>%
    mutate(w=rnorm(N, c(.3, .3)[x+1], .1))
m2 <- brm(w ~ 1 + x, data=DF2, prior = test)

plot(hypothesis(m1, "xTRUE = 0"))
plot(hypothesis(m2, "xTRUE = 0"))
head(posterior_samples(m1))
