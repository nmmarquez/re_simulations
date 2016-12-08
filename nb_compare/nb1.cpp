#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
    // Data
    DATA_ARRAY(X);
    DATA_VECTOR(y);
    
    //Parameters
    PARAMETER_VECTOR(betas);
    PARAMETER(log_alpha);
    PARAMETER_VECTOR(nu);
    
    // shortcuts
    int N = y.size(); // number of observations
    int M = betas.size(); // number of covariates
    Type alpha = exp(log_alpha);
    
    // linear model
    vector<Type> mu(N);
    Type pred;
    for(int i=0; i < N; i++){
        pred = 0.0;
        for(int b=0; b < M; b++){
            pred += (X(i, b) * betas[b]) + nu[i];
        }
        mu[i] = exp(pred);
    }
    
    // Begin calculating likelihood
    Type nll = 0.;
    
    // parameter likelihood
    for(int i=0; i < N; i++){
        nll -= dgamma(exp(nu[i]), Type(1. / alpha), alpha, true);
    }
    
    // data likelihood
    for(int i=0; i < N; i++){
        nll -= dpois(y[i], mu[i], true);
    }
    
    REPORT(alpha);
    REPORT(betas);
    REPORT(mu);
    REPORT(nu);

    return nll;
}
