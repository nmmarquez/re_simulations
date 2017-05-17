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
        Q = pow(exp(logtau), 2.) * (kappa4*M0 + Type(2.0)*kappa2*M1 + M2);
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
        DATA_VECTOR(y);
        DATA_MATRIX(cov);
        DATA_IVECTOR(geo);
        DATA_IVECTOR(temporal);
        DATA_IVECTOR(age);
        DATA_INTEGER(T);
        DATA_INTEGER(A);
        
        // SPDE objects
        DATA_SPARSE_MATRIX(M0);
        DATA_SPARSE_MATRIX(M1);
        DATA_SPARSE_MATRIX(M2);
        
        PARAMETER_VECTOR(beta);
        PARAMETER(logtau);
        PARAMETER(logkappa);
        PARAMETER_VECTOR(logitrho);
        PARAMETER(logsigma);
        PARAMETER_ARRAY(phi);
        printf("%s\n", "Parameters loaded");
        
        int N = y.size();
        int J = beta.size();
        
        vector<Type> rho = Type(1.) / (Type(1.) + exp(Type(-1.) * logitrho));
        Type sigma = exp(logsigma);
        
        SparseMatrix<Type> Q_time = ar_Q(T, rho[0], Type(1.));
        SparseMatrix<Type> Q_age = ar_Q(A, rho[1], Type(1.));
        SparseMatrix<Type> Q_loc = spde_Q(logkappa, logtau, M0, M1, M2);
        
        printf("%s\n", "Parameters built");
        
        vector<Type> yhat(N);
        
        Type nll = 0.;
        
        nll += SEPARABLE(GMRF(Q_time), SEPARABLE(GMRF(Q_age), GMRF(Q_loc)))(phi);
        
        printf("%s\n", "Random effects evaluated");
        
        for (int i = 0; i < N; i++) {
            yhat(i) = 0.;
            for (int j= 0; j < J; j++) {
                yhat(i) += beta(j) * cov(i,j);
            }
            yhat(i) += phi(geo(i), age(i), temporal(i));
            nll -= dnorm(y(i), yhat(i), sigma, true);
        }
        
        printf("%s\n", "likelihood evaluated");
        
        REPORT(phi);
        REPORT(sigma);
        REPORT(logtau);
        REPORT(logkappa);
        REPORT(rho);
        
        return(nll);
    }
