#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
using namespace density;
using Eigen::SparseMatrix;

// create a sparse precision matrix of an ar1 process
template<class Type>
SparseMatrix<Type> ar_Q(int N, Type rho, Type sigma) {
    SparseMatrix<Type> Q(N,N);
    Q.insert(0,0) = (1.) / pow(sigma, 2.);
    for (size_t n = 1; n < N; n++) {
        Q.insert(n,n) = (1. + pow(rho, 2.)) / pow(sigma, 2.);
        Q.insert(n-1,n) = (-1. * rho) / pow(sigma, 2.);
        Q.insert(n,n-1) = (-1. * rho) / pow(sigma, 2.);
    }
    Q.coeffRef(N-1,N-1) = (1.) / pow(sigma, 2.);
    return Q;
}

// create a covariance matrix of an ar1 process
template<class Type>
matrix<Type> ar_vcov(int N, Type rho, Type sigma) {
    matrix<Type> S(N,N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            S(i,j) = pow(rho, abs(i-j)) * pow(sigma, 2.) / (1-pow(rho, 2.));
        }
    }
    return S;
}

template<class Type>
Type objective_function<Type>::operator() (){
    
    DATA_VECTOR(yobs);
    DATA_INTEGER(option);
    
    PARAMETER(logsigma);
    PARAMETER(logitrho);
    
    Type rho = Type(1.) / (Type(1.) + exp(Type(-1.) * logitrho));
    Type sigma = exp(logsigma);
    Type nll = 0.;
    int N = yobs.size();
    
    // for loop method for time series
    if(option == 1){
        for(int i = 1; i < N; i++){
            nll -= dnorm(yobs(i), rho * yobs(i-1), sigma, true);
        }
    }
    
    // covariance method for time series
    if(option == 2){
        nll += MVNORM(ar_vcov(N, rho, sigma))(yobs);
    }
    
    // precision method for time series
    if(option == 3){
        nll += GMRF(ar_Q(N, rho, sigma))(yobs);
    }
    
    // built in method
    if(option == 4){
        nll += SCALE(AR1(rho), pow(pow(sigma, 2.) /(1. -  pow(rho, 2.)), .5))(yobs);
    }
    
    REPORT(rho);
    REPORT(sigma);
    
    return(nll);
}
