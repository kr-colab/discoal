/* alleleTraj.c */
/* alleleTrajectory functions */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "ranlib.h"
#include "alleleTraj.h"

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

	return(epsilon/denom);
}

double detSweepFreqGeneral(double alpha, double A) {
	double epsilon = 0.05 / alpha;
	/* For constant alpha_eff, sweep duration tau_s satisfies alpha*tau_s = -2*log(eps),
	 * so A_ts = -2*log(eps). For time-varying alpha_eff, this anchors the boundary
	 * to the user-specified alpha (interpretation: alpha is the "intended" sweep
	 * strength; the trajectory follows the time-varying alpha_eff). */
	double A_ts = -2.0 * log(epsilon);
	double denom = epsilon + ((1.0 - epsilon) * exp(A - A_ts));
	return epsilon / denom;
}


