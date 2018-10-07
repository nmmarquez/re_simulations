// Simple linear regression.
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
using namespace density;
using Eigen::SparseMatrix;

template<class Type>
SparseMatrix<Type> iid_Sigma(int N, Type sigma){
    SparseMatrix<Type> Sigma(N, N);
    for(int i = 0; i < N; i++){
        Sigma.insert(i,i) = pow(sigma, 2.);
    }
    return Sigma;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
    // Data Variables
    DATA_VECTOR(Y);
    DATA_VECTOR(x);
    DATA_FACTOR(group);
    DATA_INTEGER(constrain);
    
    // Fixed Effects
    PARAMETER(a);
    PARAMETER(b);
    PARAMETER(logSigma);
    PARAMETER(logSigmaZ);
    
    // random effects
    PARAMETER_VECTOR(z);
    
    // Predictions
    int N = z.size();
    SparseMatrix<Type> Sigma = iid_Sigma(N, exp(logSigmaZ));
    matrix<Type> A(1,N);
    vector<Type> pred(Y.size());
    vector<Type> zstar(z.size());
    matrix<Type> At;
    Type res = 0.;
    Type placeScalar = 0;
    vector<Type> placevec;
    
    
    printf("%s\n", "Constraining.");
    if(constrain == 2){
        for(int i = 0; i < N; i++){
            A(0,i) = 1.;
        }
        printf("%s\n", "Transpose A.");
        At = A.transpose();
        printf("%s\n", "Build the 1x1 matrix.");
        placeScalar = pow((A * Sigma * At)(0,0), -1);
        printf("%s\n", "Build new place vector.");
        placevec = Sigma * At;
        printf("%s\n", "Multiply by scalar.");
        placevec = placevec * placeScalar;
        printf("%s\n", "Get new 1x1 value on right most side.");
        placeScalar = (A * z)(0,0);
        printf("%s\n", "Multiply by new vector");
        placevec = placevec * placeScalar;
        printf("%s\n", "Get adjusted values.");
        zstar = z - placevec;
        printf("Final zstar size %i\n", int(zstar.size()));
        //zstar = z - (Sigma * At * pow(A * Sigma * At, -1.) * (A * z));
    }

    if(constrain == 1){
        zstar.resize(z.size() + 1);
        for(int j=0;j<z.size();j++){
            zstar[j] = z[j];
        }
        zstar[z.size()] = Type(-1.0) * sum(z);
    }

    if(constrain == 0){
        for(int j=0;j<z.size();j++){
            zstar[j] = z[j];
        }
    }
    
    for(int i=0;i<Y.size();i++){
        //printf("Making Predictions %i\n", i);
        pred[i] = a + b * x[i] + zstar[group[i]];
    }
    
    // likelihood
    for(int j=0;j<zstar.size();j++){
        if(constrain == 1){
            res -= dnorm(Type(0.0), zstar[j], exp(logSigmaZ), true);
        }
        else{
            res -= dnorm(Type(0.0), z[j], exp(logSigmaZ), true);
        }
    }
    
    printf("%s\n", "Evaluating data likelihood.");
    res -= sum(dnorm(Y, pred, exp(logSigma), true));
    return res;
}
