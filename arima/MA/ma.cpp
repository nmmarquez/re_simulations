#include <TMB.hpp>
using namespace density;
using Eigen::SparseMatrix;

// OBJECTIVE FUNCTION
template<class Type>
Type objective_function<Type>::operator() () {
    
    DATA_VECTOR(x);
    
    PARAMETER(mu);
    PARAMETER(log_sigma);
    PARAMETER(log_sigmaw);
    PARAMETER_VECTOR(q);
    PARAMETER_VECTOR(w);
    
    // name some helpful variables
    int N = x.size();
    int Nq = q.size();
    vector<Type> xhat(N);
    Type nll = 0.;
    Type sigma = exp(log_sigma);
    Type sigmaw = exp(log_sigmaw);
    
    // eval the errors
    for(int i = 0; i < N; i++){
        nll -= dnorm(Type(0.), w[i], sigmaw, true);
    }
    
    // Predict that shit 
    for(int i = Nq; i < N; i++){
        xhat[i] = mu + w[i];
        for(int j = 0; j < Nq; j++){
            xhat[i] += w[i - j - 1] * q[j];
        }
        nll -= dnorm(x[i], xhat[i], sigma, true);
    }
    
    ADREPORT(mu);
    ADREPORT(q);
    ADREPORT(sigmaw);
    REPORT(w);
    return(nll);
}
