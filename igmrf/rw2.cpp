#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
using namespace density;
using Eigen::SparseMatrix;
template<class Type>
SparseMatrix<Type> spde_Q(Type logkappa, Type logtau, SparseMatrix<Type> M0,
                          SparseMatrix<Type> M1, SparseMatrix<Type> M2) {
    SparseMatrix<Type> Q;
    Type kappa2 = exp(2. * logkappa);
    Type kappa4 = kappa2*kappa2;
    Q = pow(exp(logtau), 2.)  * (kappa4*M0 + Type(2.0)*kappa2*M1 + M2);
    return Q;
}

template<class Type>
matrix<Type> rw2_Q(int N, Type sigma) {
    matrix<Type> Q(N,N);
    Q(0,0) = (1.) / pow(sigma, 2.);
    Q(0,1) = (-2.) / pow(sigma, 2.);
    Q(1,0) = (-2.) / pow(sigma, 2.);
    Q(1,1) = (5.) / pow(sigma, 2.);
    for (int n = 2; n < (N-1); n++) {
        Q(n,n) = (6.) / pow(sigma, 2.);
        Q(n-1,n) = (-4.) / pow(sigma, 2.);
        Q(n,n-1) = (-4.) / pow(sigma, 2.);
        Q(n-2,n) = (1.) / pow(sigma, 2.);
        Q(n,n-2) = (1.) / pow(sigma, 2.);
    }
    Q(N-2,N-2) = (5.) / pow(sigma, 2.);
    Q(N-1,N-2) = (-2.) / pow(sigma, 2.);
    Q(N-2,N-1) = (-2.) / pow(sigma, 2.);
    Q(N-1,N-1) = (1.) / pow(sigma, 2.);
    return Q;
}

template<class Type>
SparseMatrix<Type> rw2_QS(int N, Type sigma) {
    Type offset = Type(.0001);
    SparseMatrix<Type> Q(N,N);
    Q.insert(0,0) = (1.) / pow(sigma, 2.) + offset;
    Q.insert(0,1) = (-2.) / pow(sigma, 2.);
    Q.insert(1,0) = (-2.) / pow(sigma, 2.);
    Q.insert(1,1) = (5.) / pow(sigma, 2.) + offset;
    for (int n = 2; n < (N-1); n++) {
        Q.insert(n,n) = (6.) / pow(sigma, 2.) + offset;
        Q.insert(n-1,n) = (-4.) / pow(sigma, 2.);
        Q.insert(n,n-1) = (-4.) / pow(sigma, 2.);
        Q.insert(n-2,n) = (1.) / pow(sigma, 2.);
        Q.insert(n,n-2) = (1.) / pow(sigma, 2.);
    }
    Q.coeffRef(N-2,N-2) = (5.) / pow(sigma, 2.) + offset;
    Q.coeffRef(N-1,N-2) = (-2.) / pow(sigma, 2.);
    Q.coeffRef(N-2,N-1) = (-2.) / pow(sigma, 2.);
    Q.insert(N-1,N-1) = (1.) / pow(sigma, 2.) + offset;
    return Q;
}

template<class Type>
Type densRW2(vector<Type>x, Type sigma) {
    Type kappa_ = pow(sigma, -2.);
    int N = x.size();
    matrix<Type> xt(1, N);
    for(int n = 1; n < (N); n++){
        xt(0, n) = x[n]; 
    }
    matrix<Type> Q = rw2_Q(N, sigma);
    matrix<Type> temp = xt * Q;
    temp = temp * x;
    Type nll = -1. * (((Type(N) - 2.) / 2.) * log(kappa_) + (-.5 * temp(0,0)));
    return nll;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
   
    DATA_VECTOR(x);
    PARAMETER(log_sigma);
    
    Type sigma = exp(log_sigma);
    
    matrix<Type> denseQ = rw2_Q(x.size(), sigma);
    
    for(int i=0;i<x.size();i++) {
        denseQ(i,i) += .00001;
    }
    
    Eigen::SparseMatrix<Type> Q = asSparseMatrix(denseQ);
    
    Type nll = GMRF(Q)(x);
    // Type nll = densRW2(x, sigma);
    
    REPORT(log_sigma);
    return nll;
}
