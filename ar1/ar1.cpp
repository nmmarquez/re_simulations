#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
using namespace density;
using Eigen::SparseMatrix;

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

template<class Type>
matrix<Type> ar_Sigma(int N, Type rho, Type sigma) {
    matrix<Type> Sigma(N,N);
    for(int i=0; i < N; i++){
        for(int j=0; j < N; j++){
            Sigma(i,j) = (sigma*sigma) / (1 - (rho*rho)) * pow(rho, abs(i-j));
        }
    }
    return Sigma;
}

template<class Type>
Type objective_function<Type>::operator() (){
    DATA_VECTOR(yobs);
    DATA_INTEGER(option);
    
    PARAMETER(beta);
    PARAMETER(log_sigma_epsilon);
    PARAMETER(log_sigma_ar);
    PARAMETER(logit_rho_ar);
    
    PARAMETER_VECTOR(zeta);
    
    printf("%s\n", "Transform parameters.");
    Type sigma_epsilon = exp(log_sigma_epsilon);
    Type sigma_ar = exp(log_sigma_ar);
    Type rho_ar = Type(1.) / (Type(1.) + exp(Type(-1.) * logit_rho_ar));
    
    int N = yobs.size();
    
    Type nll = 0.;
    
    // random effects evaluate
    
    if(option == 0){
        for(int n=1; n < N; n++){
            nll -= dnorm(zeta[n], zeta[n-1] * rho_ar, sigma_ar, true);
        }
    }
    
    if(option == 1){
        nll += GMRF(ar_Q(N, rho_ar, sigma_ar))(zeta);
    }
    
    if(option == 2){
        nll += SCALE(AR1(rho_ar), pow(sigma_ar*sigma_ar/(1-pow(rho_ar, 2)),0.5))(zeta);
    }
    
    if(option == 3){
        nll += MVNORM(ar_Sigma(N, rho_ar, sigma_ar))(zeta);
    }
    
    //fixed effects evaluate
    for(int n=1; n < N; n++){
        nll -= dnorm(yobs[n], beta + zeta[n], sigma_epsilon, true);
    }
    
    return(nll);
}
