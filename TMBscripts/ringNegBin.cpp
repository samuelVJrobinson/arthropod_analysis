#include <TMB.hpp>
/*    "Onion skin" negative binomial functional regression
 *     y ~ NegBin(mu,theta)
 *     log(mu) = (coefMat %*% coefVec) + (ringMat %*% ringDist)
 *     ringDist(d|eta,rho,lambda) = eta*exp(-(rho^2)*(d^lambda))
 */
template<class Type>
Type objective_function<Type>::operator() ()
{	
	// Data
	DATA_VECTOR(y);  //observations - counts
	DATA_MATRIX(coefMat); //design matrix for linear fixed effects
	// DATA_VECTOR(ringDist);  //vector of ring distances
	// DATA_MATRIX(ringMat); //matrix of ring compositions
  
	// Parameters
	PARAMETER_VECTOR(coefVec); //coefficients for linear fixed effects  
	PARAMETER(eta); //intercept for function
	PARAMETER(logRho); //log distance decay
	PARAMETER(lambda); //exponent of distance function
	PARAMETER(logTheta); //log theta for nbinom
	
	//Preamble
	int nCount = y.size(); //Number of observations  
	Type nll = 0; // Objective function

	//Linear component of model
	vector<Type> mu = coefMat * coefVec; 

	//Radial component of model
	vector<Type> ringCoef(ringDist.size());
	
	for(int i=0; i<ringDist.size(); i++){ //For each ring distance	
		ringCoef(i) = eta*exp(-pow(exp(logRho),2)*pow(ringDist(i),lambda));
	}

	mu = mu + (ringMat * ringCoef); //Add radial component to mu

	for(int i=0; i<nCount; i++){
	  nll = nll - dnbinom2(y,exp(mu(i)),exp(logTheta),true);
	}
	// // Reporting
	// Type theta = exp(logTheta);
	// Type rho = exp(logRho); 
	// 
	// ADREPORT(theta);
	// ADREPORT(rho);
	// ADREPORT(ringCoef);
	
	return nll;
}

