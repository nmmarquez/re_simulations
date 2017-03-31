#include <TMB.hpp>
using namespace density;
using Eigen::SparseMatrix;

// OBJECTIVE FUNCTION
template<class Type>
Type objective_function<Type>::operator() () {
    
    DATA_VECTOR(y);
    DATA_VECTOR(x);
    DATA_IVECTOR(group);
    
    PARAMETER(b0);
    PARAMETER(b1);
    PARAMETER(log_sigma);
    PARAMETER(log_re_sigma);
    PARAMETER_VECTOR(zeta);
    
    int N = y.size();
    int Z = zeta.size();
    Type sigma = exp(log_sigma);
    Type re_sigma = exp(log_re_sigma);
    vector<Type> yhat(N);
    Type nll = 0.;
    
    for (int n = 0; n < N; n++) {
        yhat[n] = b0 + b1 * x[n] + zeta[group[n]];
        nll -= dnorm(y[n], yhat[n], sigma, true);
    }
    
    for (int z = 0; z < Z; z++) {
        nll -= dnorm(zeta[z], Type(0.), re_sigma, true);
    }
    
    REPORT(b0);
    REPORT(b1);
    REPORT(sigma);
    REPORT(re_sigma);
    
    return nll;
}
