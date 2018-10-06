// Simple linear regression.
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>


template<class Type>
Type objective_function<Type>::operator() ()
{
    // Data Variables
    DATA_VECTOR(Y);
    DATA_VECTOR(x);
    DATA_FACTOR(group);
    DATA_INTEGER(constrain);
    
    // Fixed Effects
    PARAMETER(a);
    PARAMETER(b);
    PARAMETER(logSigma);
    PARAMETER(logSigmaZ);
    
    // random effects
    PARAMETER_VECTOR(z);
    
    // Predictions
    vector<Type> pred(Y.size());
    vector<Type> z_(z.size());

    if(constrain == 1){
        z_.resize(z.size() + 1);
        for(int j=0;j<z.size();j++){
            z_[j] = z[j];
        }
        z_[z.size()] = Type(-1.0) * sum(z);
    }

    if(constrain == 0){
        for(int j=0;j<z.size();j++){
            z_[j] = z[j];
        }
    }
    
    for(int i=0;i<Y.size();i++){
        pred[i] = a + b * x[i] + z_[group[i]];
    }
    
    // likelihood
    
    Type res=0;
    
    for(int j=0;j<z_.size();j++){
        res -= dnorm(Type(0.0), z_[j], exp(logSigmaZ), true);
    }
    
    res -= sum(dnorm(Y, pred, exp(logSigma), true));
    return res;
}
