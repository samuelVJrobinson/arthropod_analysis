// Ring regression
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);  //observations - counts
  DATA_MATRIX(coefMat); //design matrix for linear fixed effects

  PARAMETER_VECTOR(coefVec); //coefficients for linear fixed effects  
  PARAMETER(logTheta); //log theta for nbinom
  
  //Preamble
  int nCount = y.size(); //Number of observations 
  Type nll = 0; // Objective function
  Type theta = exp(logTheta);
  
  //Linear + radial component of model
  vector<Type> mu = coefMat * coefVec;
  
  // // Normal distribution
  // for(int i=0; i<nCount; i++){ //For each data point
  //   nll -= dnorm(y(i), mu(i), theta, true); //Increment nll
  // }
  
  // // Poisson distribution
  // for(int i=0; i<nCount; i++){ //For each data point
  //   nll -= dpois(y(i), exp(mu(i)), true); //Increment nll
  // }
  
  //Neg bin
  for(int i=0; i<nCount; i++){ //For each data point
    nll -= dnbinom2(y(i), exp(mu(i)), exp(mu(i))+(pow(exp(mu(i)),2)/theta), true); //Increment nll
  }
  
  REPORT(theta);
  ADREPORT(theta);
  
  return nll;
}
