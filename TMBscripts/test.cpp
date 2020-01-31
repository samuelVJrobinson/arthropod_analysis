#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{	
  // Data
  DATA_VECTOR(y);  //observations - counts
  DATA_MATRIX(coefMat); //design matrix for linear fixed effects
  
  // Parameters
  PARAMETER_VECTOR(coefVec); //coefficients for linear fixed effects  
  PARAMETER(logTheta); //log theta for nbinom
  
  //Preamble
  int nCount = y.size(); //Number of observations  
  Type nll = 0; // Objective function
  
  //Linear component of model
  vector<Type> mu = coefMat * coefVec;
  
  for(int i=0; i<nCount; i++){
    nll -= dnorm(y,mu(i),exp(logTheta),true);
  }
  
  return nll;
}

