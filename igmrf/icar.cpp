#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
using namespace density;
using Eigen::SparseMatrix;


template<class Type>
SparseMatrix<Type> iid_Q(int N){
    SparseMatrix<Type> Q(N, N);
    for(int i = 0; i < N; i++){
        Q.insert(i,i) = 1.;
    }
    return Q;
}

template<class Type>
SparseMatrix<Type> lcar_Q(SparseMatrix<Type> Wstar, Type rho, Type sigma){
    int N = Wstar.rows();
    
    SparseMatrix<Type> I(N, N);
    for(int i = 0; i < N; i++){
        I.insert(i,i) = 1.;
    }
    
    SparseMatrix<Type> Q = (1. / sigma) * (rho * (Wstar) + (1. - rho) * I);
    return Q;
}

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
Type objective_function<Type>::operator() (){
    
    DATA_ARRAY(yobs);
    DATA_SPARSE_MATRIX(Wstar); // pre compiled wstar matrix
    PARAMETER(log_sigmasp); // lcar effects
    PARAMETER(log_sigma);
    PARAMETER_VECTOR(zeta);
    PARAMETER(rho);

    Type spsigma = exp(log_sigmasp);
    Type sigma = exp(log_sigma);

    SparseMatrix<Type> Q_loc = lcar_Q(Wstar, rho, spsigma);
    
    // Type nll = GMRF(Q_loc)(zeta);
    Type nll = 0.;
    
    vector<Type> tmp = Wstar*zeta;
    nll -= Type(-0.5)*(zeta*tmp).sum();

    vector<Type> eta = tmp/pow(spsigma, -2.);

    for(int i=0; i<yobs.size(); i++){
        nll -= dnorm(yobs[i], eta[i], sigma, true);
    }
    
    return(nll);
}
