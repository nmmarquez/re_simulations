#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
    // Data
    DATA_ARRAY(X);
    DATA_VECTOR(y);
    
    //Parameters
    PARAMETER_VECTOR(betas);
    PARAMETER(log_theta);
    
    // shortcuts
    int N = y.size(); // number of observations
    int M = betas.size(); // number of covariates
    Type theta = exp(log_theta);
    
    // linear model
    vector<Type> mu(N);
    Type pred;
    for(int i=0; i < N; i++){
        pred = 0.0;
        for(int b=0; b < M; b++){
            pred += X(i, b) * betas[b];
        }
        mu[i] = exp(pred);
    }
    
    // Begin calculating likelihood
    Type nll = 0.;
    Type nb_var;
    
    // data likelihood
    for(int i=0; i < N; i++){
        nb_var = mu[i] + (pow(mu[i], Type(2.)) / theta);
        nll -= dnbinom2(y[i], mu[i], nb_var, true);
    }
    
    REPORT(theta);
    REPORT(betas);
    REPORT(mu);
    
    return nll;
}
