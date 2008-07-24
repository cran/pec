#include <stdlib.h> 
void survest_cox_aalen(double *hazard,
		       double *coef,
		       double *vars,
		       int *nvars,
		       int *nobs,
		       int *ntime)
{
  int t,i,z;
  /*
    this loop saves the time-varying part of the fitted hazard
  */
  for (t=0;t<*ntime;t++){ /*  at each jump time of the fit   */
    for (i=0;i<*nobs;i++){ /* and for each patient we compute the cumulative hazard  */
      for (z=0;z<*nvars;z++){ 
	hazard[i + t*(*nobs)] += coef[t + z*(*ntime)] * vars[i + z *(*nobs)];
      }
    }
  }
}
