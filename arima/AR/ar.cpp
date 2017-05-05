#include <TMB.hpp>
using namespace density;
using Eigen::SparseMatrix;


template<class Type>
vector<Type> diff(vector<Type> vec, int differences) {
    vector<Type> newvec = vec;
    vector<Type> diffvec = vec;
    for(int d = 0; d < differences; d++){
        diffvec.resize(newvec.size() - 1);
        for(int i = 0; i < diffvec.size(); i++){
            diffvec[i] = newvec[i + 1] - newvec[i];
        }
        newvec = diffvec;
    }
    return newvec;
}


// OBJECTIVE FUNCTION
template<class Type>
Type objective_function<Type>::operator() () {
    
    DATA_VECTOR(raw);
    DATA_INTEGER(d);
    //DATA_MATRIX(covs);
    
    PARAMETER(mu);
    PARAMETER(log_sigma);
    PARAMETER_VECTOR(p);
    //PARAMETER_VECTOR(b);
    
    // name some helpful variables
    vector<Type> x = diff(raw, d);
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

