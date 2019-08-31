rm(list=ls())
library(INLA)
library(sp)
library(fields)
library(geoR)
library(viridisLite)
library(tidyverse)

data('ca20')

df = data.frame(
    y = ca20$data, locx = ca20[[1]][ , 1], locy = ca20[[1]][ , 2], ca20[[3]])
spatial.scaling = 100
df$locx = (df$locx - min(df$locx))/spatial.scaling
df$locy = (df$locy - min(df$locy))/spatial.scaling
df$altitude = df$altitude - mean(df$altitude)
df$y = df$y-50

max.edge = 0.5
mesh <- inla.mesh.2d(
    loc=df[ , c('locx', 'locy')],
    offset = c(0.5, 1.5),
    max.edge=c(max.edge, max.edge*3),
    # discretization accuracy
    cutoff=max.edge/5)

A = inla.spde.make.A(mesh=mesh, loc=data.matrix(df[ , c('locx', 'locy')]))

Xcov = data.frame(intercept=1, altitude=df$altitude)
Xcov = as.matrix(Xcov)

stack <- inla.stack(tag='est',
                    # - Name (nametag) of the stack
                    # - Here: est for estimating
                    data=list(y=df$y),
                    effects=list(
                        # - The Model Components
                        s=1:mesh$n,
                        # - The "s" is means "spatial"
                        Xcov=Xcov),
                    # - The second is all fixed effects
                    A = list(A, 1)
                    # - First projector matrix is for 's'
                    # - second is for 'fixed effects'
)

plot(mesh, asp=1)
points(df[ , c('locx', 'locy')], col='red')

prior.median.sd = 1; prior.median.range = 7
spdepc = inla.spde2.pcmatern(
    mesh,
    prior.range = c(prior.median.range, .05),
    prior.sigma = c(prior.median.sd, .05), constr = T)
spde = inla.spde2.matern(
    mesh)
formula = y ~ -1 + Xcov + f(s, model=spde)
formulapc = y ~ -1 + Xcov + f(s, model=spdepc)

prior.median.gaus.sd = 5.5
family = 'gaussian'
control.family = list(hyper = list(prec = list(
    prior = "pc.prec", fixed = FALSE, param = c(prior.median.gaus.sd,0.5))))

res <- inla(formula, data=inla.stack.data(stack),
            control.predictor=list(A = inla.stack.A(stack), compute=T),
            family = family,
            control.family = control.family,
            control.inla = list(int.strategy='eb'),
            verbose=F)

respc <- inla(formulapc, data=inla.stack.data(stack),
            control.predictor=list(A = inla.stack.A(stack), compute=T),
            family = family,
            control.family = control.family,
            control.inla = list(int.strategy='eb'),
            verbose=F)


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
        
        # id 23001 gives us the spde sigma parameter untransformed
        # if we have a normal prior we backtransfrom for interpretability
        else if((prior$prior == "normal") & (id == "23001")){
            kapp <- exp(res$summary.hyperpar[i+1,"0.5quant"])
            1 / (4 * pi * exp(i.lk) ^ 2 * exp(i.lt) ^ 2)
        }
        
        # id 23002 gives us the spde kappa parameter untransformed
        # if we have a normal prior we backtransfrom for interpretability
        else if((prior$prior == "normal") & (id == "23002")){
            malt <- inla.smarginal(hh) %>%
                {inla.tmarginal(function(x) sqrt(8) / exp(x), .)} %>%
                as_tibble() %>%
                mutate(Distribution="Posterior", label=label)
            
        }
        
        else{
            # portion for any prior except pcmatern prior
            df <- bind_rows(
                m,
                INLA:::inla.get.prior.xy(
                    tolower(id[2]), id[1], all.hyper, range = range(m$x)) %>%
                    as_tibble() %>%
                    mutate(Distribution="Prior", label=label)
            )
        }
        
        filter(df, !is.na(x) & !is.na(y))
    }))
    return(hyperDF)
}

summary(res)

extract_hyperpar_prior_posteriors(res) %>%
    ggplot(aes(x=x, y=y, color=Distribution, type=Distribution)) +
    geom_line() +
    theme_classic() +
    facet_wrap(~label, scales="free")


