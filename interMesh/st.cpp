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
Type objective_function<Type>::operator() ()
{
    using namespace R_inla;
    using namespace density;
    using namespace Eigen;
    
    // Counts of observed values
    DATA_IVECTOR(yPoint);
    
    // Denoms
    DATA_IVECTOR(denomPoint);
    
    // Projections
    DATA_SPARSE_MATRIX(AprojPoint);
    
    //time periods
    DATA_INTEGER(timem);
    DATA_INTEGER(spacem);
    
    // SPDE objects
    DATA_SPARSE_MATRIX(M0);
    DATA_SPARSE_MATRIX(M1);
    DATA_SPARSE_MATRIX(M2);
    
    // Parameters
    PARAMETER(beta0);
    PARAMETER(log_tau);
    PARAMETER(log_kappa);
    PARAMETER_VECTOR(z);
    PARAMETER(logit_rho);
    
    // printf("%s\n", "Loading data complete.");
    
    int Npoint = yPoint.size();
    
    Type tau = exp(log_tau);
    Type kappa = exp(log_kappa);
    Type rho = Type(1.) / (Type(1.) + exp(Type(-1.) * logit_rho));
    
    Type nll = 0.0;
    
    // printf("%s\n", "Evaluate random effects.");
    
    SparseMatrix<Type> Qspde = spde_Q(log_kappa, log_tau, M0, M1, M2);
    //SparseMatrix<Type> Qar = ar_Q(timem, rho, Type(1.));
    
    // printf("%s\n", "Evaluating likelihood of RE latent field.");
    array<Type> zmat(spacem, timem);
    zmat << z;
    
    // Initially I tried doing this while testing out the kronecker function
    // and the speed differences were insane. A model that used seperable 
    // ran in a couple minutes while the kronecker version took like 3 hours
    // Im guessing this is because of the way that the sparseness detection
    // algorithm works with GMRF and is propegated with scale. This does not 
    // work with kronecker and thus a much large matrix must be searched for 
    // sparseness. 
    
    //nll += GMRF(kronecker(Qar, Qspde))(z);
    nll += SEPARABLE(AR1(rho), GMRF(Qspde))(zmat);
    
    vector<Type> projPoint = AprojPoint * z;
    
    for(int i=0; i<yPoint.size(); i++){
        // printf("Evaluating likelihood of Point %i\n", i);
        Type logitp = beta0 + projPoint[i];
        Type p = exp(logitp) / (Type(1) + exp(logitp));
        nll -= dbinom(Type(yPoint[i]), Type(denomPoint[i]), p, true);
    }
    
    REPORT(z);
    return nll;
}
