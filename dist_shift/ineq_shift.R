mu2015 <- 4
eps <- .4
x <- seq(mu2015 - 4*eps, 
         mu2015 + 4*eps, length=100)
hx <- dnorm(x, mean=mu2015, sd=eps)
plot(x, hx, "l", xlab="2015 Honduras Life Expectancy Forecast", ylab="Density")


plot_gam_dist <- function(sh, sc, add=""){
    x <- seq(0, 80, length.out=1000)
    hx <- dgamma(x, shape=sh, scale=sc)
    aineq <- qgamma(.95, shape=sh, scale=sc) - qgamma(.05, shape=sh, scale=sc)
    rineq <- qgamma(.95, shape=sh, scale=sc) / qgamma(.05, shape=sh, scale=sc)
    #t_ <- paste(add, "Relative:", round(rineq, 2), "Absolute", round(aineq, 2))
    cat(paste("Relative:", round(rineq, 2), "Absolute", round(aineq, 2)))
    t_ <- add
    plot (x, hx, "l", main=t_, ylab="Density", xlab="U5MR per 1000")
}

par(mfrow=c(2,2))
plot_gam_dist(24, 2, add="1990 Reference\nRelative: 2.3 Absolute: 39.18")
plot_gam_dist(48, 1/3, add="Scenario 1: Reduce Skew\nRelative: - Absolute: -")
plot_gam_dist(24, 2/3, add="Scenario 2: Reduce Var\nRelative: c Absolute: -")
plot_gam_dist(16/6, 6,
              add="Scenario 3: Increase Skew & Var\nRelative: + Absolute: c/+")
par(mfrow=c(1,1))

hx <- dgamma(x, shape=32, scale=1.5)
plot (x, hx, "l")
qgamma(.99, shape=32, scale=1.5) - qgamma(.01, shape=32, scale=1.5)
qgamma(.99, shape=32, scale=1.5) /  qgamma(.01, shape=32, scale=1.5)

x <- seq(0, 80, length.out=1000)
hx <- dgamma(x, shape=32, scale=.5)
plot (x, hx, "l")
qgamma(.99, shape=32, scale=.5) - qgamma(.01, shape=32, scale=.5)
qgamma(.99, shape=32, scale=.5) / qgamma(.01, shape=32, scale=.5)

x <- seq(0, 80, length.out=1000)
hx <- dgamma(x, shape=16/4.5, scale=4.5)
plot (x, hx, "l")
qgamma(.99, shape=sh, scale=sc) - qgamma(.01, shape=sh, scale=4.5)
qgamma(.99, shape=sh, scale=sc) / qgamma(.01, shape=sh, scale=4.5)

x <- seq(0, 80, length.out=1000)
hx <- dgamma(x, shape=16/9, scale=9)
plot (x, hx, "l")
