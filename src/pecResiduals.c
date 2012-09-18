#include <math.h>
void pecResiduals(double *pec,
		  double *resid,
		  double *Y,
		  double *D,
		  double *times,
		  double *pred,
		  double *weight,
		  double *weight_obs,
		  int *N,
		  int *NT,
		  int *cmodel,
		  int *ConstantPrediction)
{
  int s, i;
  double p, brier, gs, gi;
  
  for (s=0; s<(*NT);s++) {
    for (i=0; i<*N;i++){
      /* prediction */
      if (*ConstantPrediction==0){
	p = pred[i + s * (*N)];
      }
      else{
	p = pred[s];
      }
      /* weights */
      gs = weight[(i + s * (*N)) * (*cmodel) + s * (1-(*cmodel))];
      gi = weight_obs[i];
      if (Y[i] <= times[s])
	brier = D[i] * p * p / gi;
      else
	brier = (1-p)*(1-p) / gs;
      resid[i + s*(*N)] = brier;
      pec[s] += brier / (double) (*N);
    }
  }
}

void pecResidualsCR(double *pec,
		    double *resid,
		    double *Y,
		    double *D,
		    double *E,
		    double *times,
		    double *pred,
		    double *weight,
		    double *weight_obs,
		    int *N,
		    int *NT,
		    int *cmodel,
		    int *ConstantPrediction)
{
  int s, i;
  double p, brier, gs, gi;
  
  for (s=0; s<(*NT);s++) {
    for (i=0; i<*N;i++){
      /* prediction */
      if (*ConstantPrediction==0){
	p = pred[i + s * (*N)];
      }
      else{
	p = pred[s];
      }
      
      /* weights */
      gs = weight[(i + s * (*N)) * (*cmodel) + s * (1-(*cmodel))];
      gi = weight_obs[i];
      
      if (Y[i] <= times[s])
	brier = E[i] * D[i] * (1-p) * (1-p) / gi;
      else
	brier = p*p / gs;
      resid[i + s*(*N)] = brier;
      pec[s] += brier / (double) (*N);
    }
  }
}

