set.seed(123)
rm(list=ls())
library(splines)

time <-  seq(0.02, 10, .02)
y <- sin(time) + rnorm(length(time), sd=.2)
plot(time, y)
lines(time, lm(y ~ time)$fitted.values, col="red")

# ploynomial function
D <- 300
X <- outer(time, 1:D, "^")
plot(time, y)
lines(time, lm(y ~ X)$fitted.values, col="red")

# K splines with D degree
D <- 4
K <- 3
knots <- length(time) * (1:K)/(K+1)
X1 <- outer(time, 1:D,"^")
X2 <- outer(time, knots ,">") * outer(time, knots, "-")^D
X <- cbind(X1 ,X2)
plot(time, y)
lines(time, lm(y ~ X)$fitted.values, col="red")
