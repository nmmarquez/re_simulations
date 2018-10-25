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
    
    // Polygon index
    DATA_IVECTOR(loc);
    
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
    
    // printf("%s\n", "Loading data complete.");
    
    int Npoint = yPoint.size();
    int Npoly = yPoly.size();
    
    Type tau = exp(log_tau);
    Type kappa = exp(log_kappa);
    
    Type nll = 0.0;
    
    // printf("%s\n", "Evaluate random effects.");
    
    SparseMatrix<Type> Q = spde_Q(log_kappa, log_tau, M0, M1, M2);
    
    // printf("%s\n", "Evaluating likelihood of RE latent field.");
    nll += GMRF(Q)(z);
    
    // printf("%s\n", "Matrix Mult 1.");
    vector<Type> projPoint = AprojPoint * z;
    // printf("%s\n", "Matrix Mult 2.");
    vector<Type> projLatObs = AprojObs * z + beta0;
    // printf("%s\n", "Matrix Mult 3.");
    vector<Type> projPObs = exp(projLatObs) / (Type(1.) + exp(projLatObs));
    SparseMatrix<Type> RAprojPoly = AprojPoly.transpose();
    vector<Type> projPoly = RAprojPoly * projPObs;
    
    for(int i=0; i<yPoint.size(); i++){
        // printf("Evaluating likelihood of Point %i\n", i);
        Type logitp = beta0 + projPoint[i];
        Type p = exp(logitp) / (Type(1) + exp(logitp));
        nll -= dbinom(Type(yPoint[i]), Type(denomPoint[i]), p, true);
    }
    
    // printf("%s\n", "Evaluating likelihood of Polygons.");
    for(int i=0; i<yPoly.size(); i++){
        // printf("Evaluating likelihood of Polygon %i\n", i);
        Type p = projPoly[loc[i]];
        nll -= dbinom(Type(yPoly[i]), Type(denomPoly[i]), p, true);
    }
    
    REPORT(z);
    return nll;
}
