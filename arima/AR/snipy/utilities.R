library(data.table)
library(TMB)
dyn.load("../ar.so")

sim_ar_data <- function(N, p, mu, d=0, sigma=1){
    Np <- length(p)
    x <- rep(0, N + Np)
    mu <- ifelse(d > 0, 0, mu)
    
    for(i in 1:Np){
        x[i] <- ifelse(abs(sum(p)) < 1, mu / sum(c(1, -1 * p)), 0)
    }
    
    for(i in (Np + 1):(N+Np)){
        x[i] <- mu + sum(x[(i-1):(i-Np)] * p) + rnorm(1, sd=sigma)
    }
    
    x <- x[(Np + 1):(N + Np)]
    
    if(d > 0){
        for(j in 1:d){
            newx <- rep(0, N + 1)
            newx[2:(N+1)] <-  cumsum(x)
            x <- newx[2:(N+1)]
        }
    }
    
    return(data.table(obs=x, time=1:N))
}

run_arima_base <- function(x, M, d=0){
    results <- arima(x, c(M,d,0))
    ord <- c("intercept", paste0("ar", 1:M))[ifelse(d > 0, 2, 1):(M+1)]
    int <- results$coef["intercept"] * (1 - sum(results$coef[paste0("ar", 1:M)]))
    mu <- c(int, results$coef[1:M])[ord]
    se <- diag(results$var.coef)[ord]**.5
    terms <- c("constant", paste0("lag_", 1:M))[ifelse(d > 0, 2, 1):(M+1)]
    data.table(val=mu, se=se, type="arima base", term=terms)
}

run_arima_TMB <- function(x, M, d=0){
    model_name <- "ar"
    Params <- list(mu=0, log_sigma=0, p=rep(0, M))
    Data <- list(raw=x, d=d)
    Map <- list()
    if(d > 0){
        Map[["mu"]] <- as.factor(NA)
    }
    
    # build and optimize the objective function
    Obj <- MakeADFun(data=Data, parameters=Params, DLL=model_name, map=Map)
    Opt <- nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr)
    (sdrep <- sdreport(Obj))
    mu <- sdrep$value[ifelse(d > 0, 2, 1):(M+1)]
    se <- diag(sdrep$cov)[ifelse(d > 0, 2, 1):(M+1)]**.5
    terms <- c("constant", paste0("lag_", 1:M))[ifelse(d > 0, 2, 1):(M+1)]
    data.table(val=mu, se=se, type="arima TMB", term=terms)
}

run_models <- function(x, M, d){
    rbindlist(list(run_arima_base(x, M, d), run_arima_TMB(x, M, d)))
}

plot_models <- function(x, M, d){
    results <- run_models(x, M, d)
    ggplot(data = results, aes(x = term, y = val)) + 
        geom_point() + 
        geom_errorbar(aes(ymin = val - se*1.96,ymax = val + se*1.96)) +
        facet_wrap(~type)
}

