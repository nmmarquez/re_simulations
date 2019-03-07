rm(list=ls())
library(INLA)
n = 4000
m = 3
size = 1
sum.to.zero = FALSE

y = c()
idx = c()
jdx = c()
x = c()

## regression coofs
beta = rnorm(m)
## sum to zero?
if (sum.to.zero)
    beta = beta - mean(beta)

for(i in 1:n) {
    xx = rnorm(m, sd = 1)
    eta = beta * xx
    
    x = c(x, xx)
    ## intercept
    idx = c(idx, rep(i, m))
    ## beta
    jdx = c(jdx, 1:m)
    
    p = exp(eta)
    prob = p/sum(p)
    yy = rmultinom(1, size = size, prob = prob)
    y = c(y, yy)
}
## now we want to predict an outcome. ignore the intercept
idx = c(idx, rep(NA, m))
jdx = c(jdx, 1:m)
xx = rnorm(m, sd = 1)
x = c(x, xx)
eta.pred = beta * xx
prob.pred = exp(eta.pred)/sum(exp(eta.pred))
y = c(y, rep(NA, m))

formula = y ~ -1 +
    ## one intercept pr multinom observation
    f(idx, model="iid",
      hyper = list(
          prec = list(
              initial = log(0.000001),
              fixed = TRUE))) +
    ## the covariates: add them as the 'weights' in the second
    ## argument
    f(jdx, x, model="iid",
      ## beta's sum to zero?
      constr = sum.to.zero,
      hyper = list(
          prec = list(
              initial = log(0.001),
              fixed = TRUE)))

r = inla(formula,
         family = "poisson",
         data = data.frame(y, idx, jdx, x),
         ## need this for the prediction
         control.compute = list(config=TRUE))

print(cbind(estimate = r$summary.random$jdx$mean, truth = beta))

## do the predictons. The linear predictors we want are index n*m+1:m
## (yes, this is somewhat tricky...)
nsample = 100
xx = inla.posterior.sample(nsample, r)
## (or use target = n*m+1:m)
target = c("Predictor.6001", "Predictor.6002", "Predictor.6003")

prob = matrix(NA, nsample, 3)
for(i in 1:nsample) {
    eta = xx[[i]]$latent[target, 1]
    prob[i, ] = exp(eta) / sum(exp(eta))
}
prob.est = apply(prob,2,mean)

print(rbind(predicted = prob.est,
            truth = prob.pred,
            abs.error = abs(prob.est-prob.pred)))


