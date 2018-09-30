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
    
    // Counts of observed values
    DATA_IVECTOR(yPoint);
    DATA_IVECTOR(yPoly);
    
    // Denoms
    DATA_IVECTOR(denomPoint);
    DATA_IVECTOR(denomPoly);
    
    // Projections
    DATA_SPARSE_MATRIX(AprojPoint);
    DATA_SPARSE_MATRIX(AprojObs);
    DATA_SPARSE_MATRIX(AprojPoly);
    
    // SPDE objects
    DATA_SPARSE_MATRIX(M0);
    DATA_SPARSE_MATRIX(M1);
    DATA_SPARSE_MATRIX(M2);
    
    // Parameters
    PARAMETER(beta0);
    PARAMETER(log_tau);
    PARAMETER(log_kappa);
    PARAMETER_VECTOR(z);
    
    printf("%s\n", "Loading data complete.");
    
    Type tau = exp(log_tau);
    Type kappa = exp(log_kappa);
    
    Type nll = 0.0;
    
    printf("%s\n", "Evaluate random effects.");
    
    SparseMatrix<Type> Q = spde_Q(log_kappa, log_tau, M0, M1, M2);
    
    nll += GMRF(Q)(z);  // Negative log likelihood
    
    printf("%s\n", "Project Points.");
    vector<Type> projPoint = AprojPoint * z;
    printf("%s\n", "Project observed.");
    vector<Type> projObs = AprojObs * z;
    printf("%s\n", "transpose matrix.");
    SparseMatrix<Type> RAprojPoly = AprojPoly.transpose();
    printf("%s\n", "Project polygon values.");
    vector<Type> projPoly = RAprojPoly * projObs;
    
    for(int i=0; i<yPoint.size(); i++){    
        Type logitp = beta0 + projPoint[i];
        Type p = exp(logitp) / (Type(1) + exp(logitp));
        nll -= dbinom(Type(yPoint[i]), Type(denomPoint[i]), p, true);
    }
    
    REPORT(z);
    return nll;
}
