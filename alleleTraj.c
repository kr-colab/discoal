/* alleleTraj.c */
/* alleleTrajectory functions */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "ranlib.h"
#include "alleleTraj.h"

//coth- hyperbolic cotangent
double coth(double x){
	double top, bottom;
	top = exp(2.0 * x) + 1.0;
	bottom = exp(2.0 * x) - 1.0;
	return(top / bottom);
}

/*detSweepFreq-- returns the frequency of the beneficial
mutation at time (0 <= t <= ts) where ts is an arbitrary
	constant. Assumes sweeps are deterministic according to
	Stephan et al. 1992 */
double detSweepFreq(double t, double alpha){
	double epsilon, ts;
	double denom;
	epsilon = 0.05/alpha;
	ts = -2.0 * log(epsilon) / alpha;
	denom = epsilon + ((1.0 - epsilon)*exp(alpha * (t - ts)));
//	printf("alpha: %e t: %e denom: %e epsilon: %e return: %e \n",alpha,t,denom,epsilon,epsilon/denom);

	return(epsilon/denom);
}

/* neutralStochastic-- returns the frequency of a neutral allele
which is sweeping through the population. This is the jump process
corresponding to the condition diffusion towards loss (i.e. backwards).
takes the dt, small time increment, and the current allele freq (i.e. 1 for fixation) */
double neutralStochastic(double dt, double currentFreq){
	// Precompute common expressions
	double drift = -currentFreq * dt;
	double p_q = currentFreq * (1.0 - currentFreq);
	
	// Early return for boundary cases to avoid sqrt of negative
	if (p_q <= 0.0) {
		return currentFreq + drift;
	}
	
	double diffusion = sqrt(p_q * dt);
	
	// Use arithmetic instead of branching: sign = ±1
	double sign = 2.0 * ranf() - 1.0;
	
	return currentFreq + drift + sign * diffusion;
}

/* Note: Optimized versions are now in alleleTraj.h as inline functions */

/* genicSelectionStochastic-- returns the frequency of a selected allele
which is sweeping through the population. This is the jump process
corresponding to the condition diffusion towards loss (i.e. backwards).
takes the dt, small time increment, the current allele freq (i.e. 1 for fixation),
and alpha the selection coefficient */
double genicSelectionStochastic(double dt, double currentFreq, double alpha){
	double newFreq;

	double ux =  (-0.5 * alpha * currentFreq * (1.0-currentFreq))* fast_coth(0.5 * alpha * currentFreq * (1.0-currentFreq));
	if (ranf() < 0.5){
		newFreq = (ux * dt) + sqrt(currentFreq * (1.0 - currentFreq) * dt);
	}
	else{
		newFreq = (ux * dt) - sqrt(currentFreq * (1.0 - currentFreq) * dt);
	}
	return(newFreq);
}

/* genicSelectionStochastic-- returns the frequency of a selected allele
which is sweeping through the population. This is the jump process
corresponding to the condition diffusion towards fixation (i.e. forwards).
takes the dt, small time increment, the current allele freq (i.e. 1 for fixation),
and alpha the selection coefficient */
double genicSelectionStochasticForwards(double dt, double currentFreq, double alpha){

	double ux = (alpha*currentFreq*(1.-currentFreq))/fast_tanh(alpha*currentFreq);
	if (ranf() < 0.5){
	 	currentFreq += (ux * dt) + sqrt(currentFreq * (1.0 - currentFreq) * dt);
	}
	else{
		currentFreq +=  (ux * dt) - sqrt(currentFreq * (1.0 - currentFreq) * dt);
	}
	return(currentFreq);
}

/*variablePopnSizeTraj -- uses the method of Takahata & Kimura to generate
conditional trajectory in variable population size */
double variablePopnSizeTraj(double dt, double currentFreq, double alpha, double h, double f){
	double mean, var, p;
	
	if(alpha > 0.0)
		mean = -1.*alpha*currentFreq*(1.-currentFreq)*(currentFreq+h*(1.-2.*currentFreq));
	else
		mean = -1.0*currentFreq;
		
	var = currentFreq*(1.-currentFreq)*2./f;
	p = (mean*dt) + ((2.*ranf()-1.)*sqrt(3.*var) *dt);
	//printf("e:%f v:%f p:%f pCur: %f out:%f\n",mean,var,p,currentFreq,p+currentFreq);
	if((currentFreq + p) > 1.0){
		return(1.0);
	}
	else{
		return(currentFreq + p );
	}
}




