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
Type objective_function<Type>::operator() ()
{
    using namespace R_inla;
    using namespace density;
    using namespace Eigen;
    
    DATA_IVECTOR(y);
    DATA_VECTOR(x);    
    DATA_IVECTOR(geo);
    // SPDE objects
    DATA_SPARSE_MATRIX(M0);
    DATA_SPARSE_MATRIX(M1);
    DATA_SPARSE_MATRIX(M2);
    DATA_IVECTOR(denom);
    
    PARAMETER(beta0);
    PARAMETER(beta1);
    PARAMETER(log_tau);
    PARAMETER(log_kappa);
    PARAMETER_VECTOR(z);
    
    Type tau = exp(log_tau);
    Type kappa = exp(log_kappa);
    
    Type nll = 0.0;
    
    SparseMatrix<Type> Q = spde_Q(log_kappa, log_tau, M0, M1, M2);;
    
    nll += GMRF(Q)(z);  // Negative log likelihood
    
    for(int i=0; i<y.size(); i++){    
        Type logitp = beta0 + beta1 * x[i] + z[geo[i]];
        Type p = exp(logitp) / (Type(1) + exp(logitp));
        nll -= dbinom(Type(y[i]), Type(denom[i]), p, true);
    }
    
    REPORT(z);
    return nll;
}
