// Simple linear regression.
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);  //observations - counts
  DATA_MATRIX(coefMat); //design matrix for linear fixed effects
  DATA_VECTOR(ringDist);  //vector of ring distances
  DATA_MATRIX(ringMat); //matrix of ring composition
  
  PARAMETER_VECTOR(coefVec); //coefficients for linear fixed effects  
  PARAMETER(eta); //intercept for function
  PARAMETER(logRho); //log distance decay
  PARAMETER(lambda); //exponent of distance function
  PARAMETER(logTheta); //log theta for nbinom
  
  //Preamble
  int nCount = y.size(); //Number of observations 
  int nRing = ringDist.size(); //Number of rings
  Type nll = 0; // Objective function
  
  //Radial component of model
  vector<Type> ringCoef(nRing);
  
  for(int i=0; i<nRing; i++){ //For each ring distance	
    ringCoef(i) = eta*exp(-pow(exp(logRho),2)*pow(ringDist(i),lambda));
  }
  
  //Linear + radial component of model
  vector<Type> mu = (coefMat * coefVec) + (ringMat * ringCoef);
  
  for(int i=0; i<nCount; i++){ //For each data point	
    nll -= dnbinom2(y(i), exp(mu(i)), exp(logTheta), true); //Increment nll
  }
  return nll;
}
