rm(list=ls())
library(INLA)
library(sp)
library(fields)
library(geoR)
library(viridisLite)
library(tidyverse)

time <- 1:6
loc <- as.matrix(expand.grid(seq(0, 3, 0.2), seq(0, 3, 0.2)))
locdist <- as.matrix(dist(loc))

## Space covariance -- sigma = 1, range = 0.3
## Matern covariance # to be compared with inla: smoothness(nu) = alpha - d/2, d = 2 here and default in inla alpha = 2.
Vs <- Matern(d =locdist, range = 0.3, smoothness = 1, phi = 1) 


## Time covariance -- r = 0.8
Vt <- diag(6) 
Vt <- 1* 0.8^abs(row(Vt)-col(Vt))

## Cross covariance
Vc <- kronecker(Vs, Vt)

## simulate the data
set.seed(10)
xx <- crossprod(chol(Vc), rep(rnorm(nrow(loc) * 6)))

## Create the time spatial data frame
simdf <- data.frame(x = rep(loc[,1], each = 6), y = rep(loc[,2], each = 6), 
                    dat= xx, datn = xx + 0.5*rnorm(length(xx)), time = rep(1:6,  nrow(loc)) )

## plot
lattice::levelplot(datn~ x + y|time, data = simdf)

knots = seq(1, 6, length = 6)
mesh1 = inla.mesh.1d(loc = knots, degree = 2, boundary = "free")

## generate space mesh 
locs <- unique(simdf[, c("x", "y")])
mesh2 <- inla.mesh.2d(loc = locs, offset = c(0.1, 0.4), max.edge = 0.3 )

## prior parameters
range0 <- 0.1
sigma0 <- 1
lkappa0 <- log(8)/2 - log(range0)
ltau0 <- 0.5*log(1/(4*pi)) - log(sigma0) - lkappa0

## build the spatial spde
spde <- inla.spde2.matern(mesh2, B.tau = matrix(c(ltau0, -1, 1),1,3),
                          B.kappa = matrix(c(lkappa0, 0, -1), 1,3), 
                          theta.prior.mean = c(0,0), theta.prior.prec = c(0.1, 1))


## build the space time indices
index = inla.spde.make.index("space", n.spde = spde$n.spde, n.group = mesh1$m)

## Link data and process 
A <- inla.spde.make.A(mesh2, loc = cbind(simdf$x, simdf$y), 
                      group = simdf$time, group.mesh = mesh1)

miniA <- inla.spde.make.A(mesh2, loc = cbind(simdf$x, simdf$y))



stack = inla.stack(
  data = list(y = simdf$datn), 
  A = list(A, 1),
  effects = list(
    space = index,
    data.frame(INT=rep(1, nrow(simdf)))))

simp = "expression:\n  logdens = -log_precision/2;\n  return(logdens)"
simpprior = list(rho=list(prior=simp))
    
## the model
formula = y ~ -1 + 
    f(space, model = spde, group = space.group,
      control.group = list(
          model = "ar1",
          hyper = list(theta1 = list(prior="pc.cor0", param=c(0.9, 0.9)) )))

res <- inla(formula, data=inla.stack.data(stack), 
            control.predictor=list(compute=TRUE, A=inla.stack.A(stack)))

`inla.all.hyper.postprocess` = function(all.hyper)
{
    ## postprocess all.hyper, by converting and replacing prior = 'mvnorm' into its p
    ## marginals. this is for the spde-models
    
    len.n = function(param, max.dim = 10000)
    {
        len = function(n) {
            return (n + n^2)
        }
        
        len.target = length(param)
        for(n in 1:max.dim) {
            if (len(n) == len.target) {
                return (n)
            }
        }
        stop(paste("length(param) is wrong:", len.target))
    }
    
    get.mvnorm.marginals = function(param)
    {
        n = len.n(param)
        mu = param[1:n]
        Q = matrix(param[-(1:n)], n, n)
        Sigma = solve((Q + t(Q))/2.)
        return (list(mean = mu, prec = 1/diag(Sigma)))
    }
    
    for (i in seq_along(all.hyper$random)) {
        for(j in seq_along(all.hyper$random[[i]]$hyper)) {
            
            if (all.hyper$random[[i]]$hyper[[j]]$prior == "mvnorm") {
                ## replace this one, and the p-following ones, with its marginals
                m = get.mvnorm.marginals(all.hyper$random[[i]]$hyper[[j]]$param)
                for(k in 1:length(m$mean)) {
                    kk = j + k - 1
                    all.hyper$random[[i]]$hyper[[kk]]$prior = "normal"
                    all.hyper$random[[i]]$hyper[[kk]]$param = c(m$mean[k], m$prec[k])
                }
            }
        }
    }
    
    return (all.hyper)
}

extract_hyperpar_prior_posteriors <- function(res){
    all.hyper <- inla.all.hyper.postprocess(res$all.hyper)
    hyper <- res$marginals.hyperpar
    nhyper <- length(hyper)
    # loop through all hyperpriors
    hyperDF <- bind_rows(lapply(1:nhyper, function(i){
        hh <- hyper[[i]]
        label <- INLA:::inla.nameunfix(names(hyper)[i])
        # Extract the posterior valuse into m
        m <- as_tibble(inla.smarginal(hh)) %>%
            mutate(Distribution="Posterior", label=label)
        id <- unlist(strsplit(attr(hyper[[i]], "hyperid"), "\\|"))
        prior <- INLA:::inla.extract.prior(tolower(id[2]), id[1], all.hyper)
        print(prior)
        # Extract the prior values, if else statements required for pcmatern
        if(prior$prior == "none"){
            df <- m
        }
        
        else if(prior$prior == "pcmatern"){
            # extract lamda1, lamda2 and d from inla
            l1 <- INLA:::inla.extract.prior(id[2], id[1], all.hyper)$param[1]
            l2 <- INLA:::inla.extract.prior(id[2], id[1], all.hyper)$param[2]
            d <- INLA:::inla.extract.prior(id[2], id[1], all.hyper)$param[3]
            
            # caalculate pc priors as done in the Fuglstad 2017 paper
            calc_pc_prior_ranges <- function(ranges){
                d/2 * l1 * ranges^(-d/2 - 1) * exp(-l1 * ranges^(-d/2))
            }
            
            calc_pc_prior_sigmas <- function(sigmas){
                l2 * exp(- l2 * sigmas)
            }
            
            m2 <- as_tibble(inla.smarginal(hyper[[i+1]]))
            offset_ <- m2$x[2] - m2$x[1]
            df <- bind_rows(
                m,
                m %>%
                    mutate(y=calc_pc_prior_ranges(x)) %>%
                    mutate(Distribution="Prior", label=label),
                tibble(x=seq(offset_, max(m2$x), by=offset_)) %>%
                    mutate(y=calc_pc_prior_sigmas(x)) %>%
                    mutate(
                        Distribution="Prior",
                        label=gsub("Range", "Stdev", label))
            )
        }
        
        else{
            # portion for any prior except pcmatern prior
            
            r_ <- range(m$x)
            df <- bind_rows(
                m,
                INLA:::inla.get.prior.xy(
                    tolower(id[2]), id[1], all.hyper, range = r_) %>%
                    as_tibble() %>%
                    mutate(Distribution="Prior", label=label)
            )
        }
        
        filter(df, !is.na(x) & !is.na(y))
    }))
    return(hyperDF)
}

summary(res)

# gamma priors on log precision are super flat!!!
# I dont know why this is the default versus some exponential distribution.
extract_hyperpar_prior_posteriors(res) %>%
    ggplot(aes(x=x, y=y, color=Distribution, type=Distribution)) +
    geom_line() +
    theme_classic() +
    facet_wrap(~label, scales="free")
