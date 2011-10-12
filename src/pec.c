#include <math.h>
#include <R.h>
/* survival probabilities */
void pec(double *pec,
	 double *Y,
	 double *D,
	 double *times,
	 double *pred,
	 double *weight,
	 double *weight_obs,
	 int *N,
	 int *NT,
	 int *cmodel,
	 int *pmodel)
{
  int s, i;
  double p, brier, gs, gi;
  
  for (s=0; s<(*NT);s++) {
    for (i=0; i<*N;i++){
      
      /* prediction */
      if (*pmodel==1)
	p = pred[i + s * (*N)];
      else
	p = pred[s];
      
      /* weights */
      gs = weight[(i + s * (*N)) * (*cmodel) + s * (1-(*cmodel))];
      gi = weight_obs[i];
      
      if (Y[i] <= times[s])
	brier = D[i] * p * p / gi;
      else
	brier = (1-p)*(1-p) / gs;
      pec[s] += brier / (double) (*N);
    }
  }
}

/* event probabilities - competing risks */

void pecCR(double *pec,
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
	   int *pmodel)
{
  int s, i;
  double p, brier, gs, gi;
  
  for (s=0; s<(*NT);s++) {
    for (i=0; i<*N;i++){
      
      /* prediction */
      if (*pmodel==1)
	p = pred[i + s * (*N)];
      else
	p = pred[s];
      
      /* weights */
      gs = weight[(i + s * (*N)) * (*cmodel) + s * (1-(*cmodel))];
      gi = weight_obs[i];
      
      if (Y[i] <= times[s])
	/* brier =  (D[i] * (E[i]-p) * (E[i]-p)) / gi; */
	if (E[i]==1)
	  brier =  (D[i] * (1-p) * (1-p)) / gi;
	else
	  brier = (D[i] * p * p) / gi;
      else
	brier = p*p / gs;

      pec[s] += brier / (double) (*N);
      /*       Rprintf("i=%d\tY[i]=%1.2f\ttimes[s]=%1.2f\tE[i]=%1.2f\tD[i]=%1.2f\tp=%1.2f\tbrier=%1.2f\tpec[s]=%1.2f\tgi=%1.2f\tgs=%1.2f\n",i,Y[i],times[s],E[i],D[i],p,brier,pec[s],gi,gs);    */
    }
  }
}


void pec_uncens(double *pec,
		double *Y,
		double *times,
		double *pred,
		int *N,
		int *NT,
		int *pmodel,
		int *survP)
{
  int s, i;
  double p, brier;
  
  for (s=0; s<(*NT);s++) {
    for (i=0; i<*N;i++){
      /* prediction */
      if (*pmodel==1)
	p = pred[i + s * (*N)];
      else
	p = pred[s];

      if (*survP==1) 
	if (Y[i] <= times[s])
	  brier = p * p;
	else
	  brier = (1-p)*(1-p);
      else
	if (Y[i] > times[s])
	  brier = p * p;
	else
	  brier = (1-p)*(1-p);
      
      pec[s] += brier / (double) *N;
    }
  }
}

void pec_noinf(double *pec,
	       double *Y,
	       double *D,
	       double *times,
	       double *pred,
	       double *weight,
	       double *weight_obs,
	       int *N,
	       int *NT,
	       int *cmodel,
	       int *pmodel)
{
  int s, i, j;
  double p, brier, gs, gi;

  for (s=0; s<*NT;s++) {
    for (j=0; j<*N; j++){
      
      /* prediction */
      p = pred[(j + s * (*N)) * (*pmodel) + s * (1-(*pmodel))];
      
      for (i=0; i<(*N); i++){
	/* weights */
	gs = weight[(i + s * (*N)) * (*cmodel) + s * (1-(*cmodel))];
	gi = weight_obs[i];
	if (Y[i] <= times[s])
	  brier = D[i] * p * p / gi;
	else
	  brier = (1-p)*(1-p) / gs;
	pec[s] += brier / (double) ((*N) * (*N));
      }
    }
  }
}
  
void pec_noinfCR(double *pec,
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
		 int *pmodel)
{
  int s, i, j;
  double p, brier, gs, gi;

  for (s=0; s<*NT;s++) {
    for (j=0; j<*N; j++){
      
      /* prediction */
      p = pred[(j + s * (*N)) * (*pmodel) + s * (1-(*pmodel))];
      
      for (i=0; i<(*N); i++){
	/* weights */
	gs = weight[(i + s * (*N)) * (*cmodel) + s * (1-(*cmodel))];
	gi = weight_obs[i];
	if (Y[i] <= times[s])
	  brier = E[i] * D[i] * (1-p) * (1-p) / gi;
	else
	  brier = p*p / gs;
	pec[s] += brier / (double) ((*N) * (*N));
      }
    }
  }
}

void pec_cmprsk(double *pec,
		double *Y,
		double *D,
		double *times,
		double *pred,
		double *weight,
		double *weight_obs,
		int *N,
		int *NT,
		int *cmodel,
		int *pmodel)
{
  int s, i;
  double p, brier, gs, gi;
  
  for (s=0; s<(*NT);s++) {
    
    for (i=0; i<(*N);i++){
      
      /* prediction */
      if (*pmodel==1)
	p = pred[i + s * (*N)];
      else
	p = pred[s];
      
      /* weights */
      gs = weight[(i + s * (*N)) * (*cmodel) + s * (1-(*cmodel))];
      gi = weight_obs[i];
      
      if (Y[i] <= times[s] && D[i]==1){
	brier = (p * p) + (1 - 2 * p)/gi;}
      else 
	brier = (p * p);
      
      pec[s] += brier / (double) (*N);
    }
  }
}
