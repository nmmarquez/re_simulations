#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
using namespace density;
using Eigen::SparseMatrix;

template<class Type>
Type transRho(Type logrho_) {
    Type rho_ = ((1. / (1 + exp(-logrho_))) - .5)*2;
    
    return rho_;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
    // Data
    DATA_IVECTOR(group);
    DATA_ARRAY(covs);
    DATA_IVECTOR(cluster);
    DATA_INTEGER(C); // number of clusters
    
    //Parameters
    PARAMETER_ARRAY(betas); // M x (G-1)
    PARAMETER_ARRAY(zetas); // C x (G -1)
    PARMETER_VECTOR(sigmas); // 2 one for each alternative group
    PARAMETER(logrho); // unbounded correlation coefficient
    
    // Parameter transforms
    Type rho = transRho(logrho);
    
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
            pred = zetas[cluster[i], g-1];
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
    vector<Type> ind_prob(G);
    vector<Type> obs_val(G);
    
    for(int i=0; i < N; i++){
        for(int g=0; g < G; g++){
            ind_prob[g] = prob(i,g);
            if (group[i] == g){
                obs_val[g] = Type(1);
            }
            else{
                obs_val[g] = Type(0);  
            }
        }
        nll -= dmultinom(obs_val, ind_prob, true);
    }
    
    // do the random effects calcs
    matrix<Type> Sigma(2,2);
    Sigma(0,0) = pow(sigmas[0], 2.);
    Sigma(1,1) = pow(sigmas[0], 2.);
    Sigma(1,0) = sigmas[1] * sigmas[2] * rho;
    Sigma(0,1) = sigmas[1] * sigmas[2] * rho;
    
    for(int c=0; c<C; c++){
        nll += MVNORM(S)(zetas(c,));
    }
    
    REPORT(prob);
    REPORT(betas);
    
    return nll;
}
