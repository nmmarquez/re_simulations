library(data.table)
library(TMB)
dyn.load("../ar.so")

sim_ar_data <- function(N, p, mu, sigma=1){
    Np <- length(p)
    x <- rep(0, N)
    
    for(i in 1:Np){
        if(abs(sum(p)) < 1){
            x[i] <- mu / sum(c(1, -1 * p))
        }
        else{
            x[i] <- mu * i
        }
    }
    
    for(i in (Np + 1):N){
        x[i] <- mu + sum(x[(i-1):(i-Np)] * p) + rnorm(1, sd=sigma)
    }
    return(data.table(obs=x, time=1:N))
}

run_arima_base <- function(x, M){
    results <- arima(x, c(Np,0,0))
    mu <- c(results$coef[M+1], results$coef[1:M])
    se <- c(diag(results$var.coef)[M+1], diag(results$var.coef)[1:M])**.5
    terms <- c("constant", paste0("lag_", 1:M))
    data.table(val=mu, se=se, type="arima base", term=terms)
}

run_arima_TMB <- function(x, M){
    model_name <- "ar"
    Params <- list(mu=0, log_sigma=0, p=rep(0, M))
    Data <- list(x=x)
    
    # build and optimize the objective function
    Obj <- MakeADFun(data=Data, parameters=Params, DLL=model_name)
    system.time(Opt <- nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr))
    sdrep <- sdreport(Obj)
    mu <- sdrep$value
    se <- diag(sdrep$cov)**.5
    terms <- c("constant", paste0("lag_", 1:M))
    data.table(val=mu, se=se, type="arima TMB", term=terms)
}

run_models <- function(x, M){
    rbindlist(list(run_arima_base(x, M), run_arima_TMB(x, M)))
}

plot_models <- function(x, M){
    results <- run_models(x, M)
    ggplot(data = results, aes(x = term, y = val)) + 
        geom_point() + 
        geom_errorbar(aes(ymin = val - se*1.96,ymax = val + se*1.96)) +
        facet_wrap(~type)
}

