#include <TMB.hpp>
using namespace density;
using Eigen::SparseMatrix;

// OBJECTIVE FUNCTION
template<class Type>
Type objective_function<Type>::operator() () {
    
    DATA_VECTOR(x);
    
    PARAMETER(mu);
    PARAMETER(log_sigma);
    PARAMETER_VECTOR(p);
    
    // name some helpful variables
    int N = x.size();
    int Np = p.size();
    vector<Type> xhat(N);
    Type nll = 0.;
    Type sigma = exp(log_sigma);
    
    // Predict that shit 
    for(int i = 0; i < N; i++){
        if(i < Np){
            xhat[i] = mu / (1 + (-1. * p.sum()));
        }
        else{
            xhat[i] = mu;
            for(int j = 0; j < Np; j++){
                xhat[i] += x[i - j - 1] * p[j];
            }
            nll -= dnorm(x[i], xhat[i], sigma, true);
        }
    }
    
    ADREPORT(mu);
    ADREPORT(p);
    REPORT(xhat);
    REPORT(Np);
    return(nll);
}

