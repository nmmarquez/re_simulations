#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
    // Data
    DATA_IVECTOR(group);
    DATA_ARRAY(covs);
    
    //Parameters
    PARAMETER_ARRAY(betas);
    
    // shortcuts
    int N = group.size(); // number of observations
    int M = betas.dim(0); // number of betas
    int G = betas.dim(1) + 1; // number of groups 
    
    // Predictions
    // get the odds ratio from the betas
    array<Type> odds_ratio(N, G-1);
    Type pred;
    for(int g=1; g < G; g++){
        for(int i=0; i < N; i++){
            pred = 0.0;
            for(int b=0; b < M; b++){
                pred += covs(i, b) * betas(b, g-1);
            }
            odds_ratio(i, g-1) = exp(pred);
        }
    }
    
    // back calculate the individual probabilities
    Type p;
    array<Type> prob(N, G);
    for(int i=0; i < N; i++){
        p = 1.;
        for(int g=1; g < G; g++){
            p += odds_ratio(i, g-1);
        }
        prob(i,0) = pow(p, -1.);
        for(int g=1; g < G; g++){
            prob(i,g) = odds_ratio(i, g-1) * prob(i,0);
        }
    }
    
    // calculate the likelihood
    Type nll = 0.;
    
    for(int i=0; i < N; i++){
        for(int g=1; g < G; g++){
            if (group[i] == g){
                nll -= dbinom(Type(1.), Type(1.), prob(i,g), true);
            }
            else{
                nll -= dbinom(Type(0.), Type(1.), prob(i,g), true);  
            }
        }
    }
    
    REPORT(prob);
    REPORT(betas);
    
    return nll;
}
