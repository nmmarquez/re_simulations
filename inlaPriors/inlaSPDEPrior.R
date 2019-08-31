rm(list=ls())
library(gridExtra)
library(lattice)
library(INLA)
library(splancs)
library(fields)
library(tidyverse)
data(PRprec)
data(PRborder)

Y <- rowMeans(PRprec[,3+1:31])
ind <- !is.na(Y)
Y <- Y[ind]
coords <- as.matrix(PRprec[ind,1:2])
alt <- PRprec$Altitude[ind]

prdomain <- inla.nonconvex.hull(coords, -0.03, -0.05,resolution=c(100,100))
prmesh <- inla.mesh.2d(boundary=prdomain, max.edge=c(.45,1),cutoff=0.2)
plot(prmesh, asp=1, main=""); lines(PRborder, col=3)
points(coords[,1], coords[,2], pch=19, cex=.5, col="red")
A <- inla.spde.make.A(prmesh, loc=coords)
seaDist <- apply(spDists(coords, PRborder[1034:1078,],longlat=TRUE), 1, min)

spde <- inla.spde2.matern(prmesh, alpha=2)
mesh.index <- inla.spde.make.index(name="field",n.spde=spde$n.spde)
stk.dat <- inla.stack(
    data=list(y=Y), A=list(A,1), tag="est",
    effects=list(
        c(mesh.index,list(Intercept=1)),
        list(
            long=inla.group(coords[,1]),
            lat=inla.group(coords[,2]),
            seaDist=inla.group(seaDist))))

spdepc <- inla.spde2.pcmatern(
    prmesh, alpha=2,
    prior.range = c(3, .05),
    prior.sigma = c(1, .05))
mesh.index.pc <- inla.spde.make.index(name="field",n.spde=spdepc$n.spde)
stk.dat.pc <- inla.stack(
    data=list(y=Y), A=list(A,1), tag="est",
    effects=list(
        c(mesh.index.pc,list(Intercept=1)),
        list(
            long=inla.group(coords[,1]),
            lat=inla.group(coords[,2]),
            seaDist=inla.group(seaDist))))


f.s <- y ~ -1 + Intercept + f(field, model=spde) +
    f(seaDist, model="rw1")
f.s.pc <- y ~ -1 + Intercept + f(field, model=spdepc) +
    f(seaDist, model="rw1")

f.s <- y ~ -1 + Intercept + f(field, model=spde) +
    f(seaDist, model="rw1",
      hyper = list(prec = list(prior="pc.prec", fixed=F, param=c(5,0.5))))
f.s.pc <- y ~ -1 + Intercept + f(field, model=spdepc) +
    f(seaDist, model="rw1",
      hyper = list(prec = list(prior="pc.prec", fixed=F, param=c(5,0.5))))

res <- inla(
    f.s, family="Gamma",data=inla.stack.data(stk.dat), verbose=TRUE,
    control.predictor=list(A=inla.stack.A(stk.dat),compute=TRUE))
respc <- inla(
    f.s.pc, family="Gamma",data=inla.stack.data(stk.dat.pc), verbose=TRUE,
    control.predictor=list(A=inla.stack.A(stk.dat.pc),compute=TRUE))

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
        
        # id 23001 gives us the spde sigma parameter untransformed if we
        # also have a normal prior we backtransfrom for interpretability
        else if((prior$prior == "normal") & (id[1] == "23001")){
            # TODO: only using the poster median kappa for this calc right now
            kapp <- exp(res$summary.hyperpar[i+1,"0.5quant"])
            fx <- function(x) (1 / (4 * pi * kapp^2 * exp(x) ^ 2))^.5
            label <- gsub("Theta1", "Stdev", label)
            malt <- inla.smarginal(hh) %>%
                {inla.tmarginal(fx, .)} %>%
                as_tibble() %>%
                mutate(Distribution="Posterior", label=label)
            r_ <- range(m$x)
            df <- bind_rows(
                malt,
                INLA:::inla.get.prior.xy(
                    tolower(id[2]), id[1], all.hyper, range = r_) %>%
                    {inla.tmarginal(fx, .)} %>%
                    as_tibble() %>%
                    mutate(Distribution="Prior", label=label) %>%
                    filter(x < max(malt$x))
            )
        }
        
        # id 23002 gives us the spde kappa parameter untransformed if we
        # also have a normal prior we backtransfrom for interpretability
        else if((prior$prior == "normal") & (id[1] == "23002")){
            label <- gsub("Theta2", "Range", label)
            malt <- inla.smarginal(hh) %>%
                {inla.tmarginal(function(x) sqrt(8) / exp(x), .)} %>%
                as_tibble() %>%
                mutate(Distribution="Posterior", label=label)
            r_ <- range(m$x)
            df <- bind_rows(
                malt,
                INLA:::inla.get.prior.xy(
                    tolower(id[2]), id[1], all.hyper, range = r_) %>%
                    {inla.tmarginal(function(x) sqrt(8) / exp(x), .)} %>%
                    as_tibble() %>%
                    mutate(Distribution="Prior", label=label) %>%
                    filter(x < max(malt$x))
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

bind_rows(
    extract_hyperpar_prior_posteriors(res) %>%
        mutate(prior="Standard Priors"),
    extract_hyperpar_prior_posteriors(respc) %>%
        mutate(prior="PC Priors")) %>%
    group_by(prior, label) %>%
    mutate(y=y/max(y)) %>%
    ggplot(aes(x=x, y=y, color=Distribution, type=Distribution)) +
    geom_line() +
    theme_classic() +
    facet_grid(prior~label, scales="free") +
    theme(
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_text(angle=90, hjust=1))


    


