#include <R.h>
void auc(double *AUC,
	 double *conc,
	 double *pairs,
	 int *tindex,
	 double *Y,
	 int *status,
	 double *times,
	 double *weight_i,/* G(T|X) */
	 double *weight, /* G(s|X) */
	 double *pred,
	 int *N,
	 int *NT,
	 int *tiedpredIn,
	 int *cens_model){
  int i,j,s;
  double wi, wj, ww;
  for (s=0; s<(*NT);s++) {
    conc[s]=0;
    pairs[s]=0;
    for (i=0;i<(*N);i++){
      /*
	for usuable pairs the smaller time must be uncensored
      */
      if (Y[i]<=times[s] && status[i]==1){
	for (j=tindex[s];j<*N;j++){
	  /*
	    censoring survival weights:
	    G(T_i-|X_i) for i
	    G(s|X_j)  for j
	  */
	  wi = weight_i[i];
	  wj = weight[(j + s * (*N)) * (*cens_model) + s * (1-(*cens_model))];
	  ww = (wi * wj);
	  /*
	    pair unusuable if any weight==0
	  */
	  if (wj>0 && wi>0){
	    if (pred[i + s * (*N)] == pred[j + s * (*N)]) {
	      /*
		call pair unusuable if same predictions
	      */
	      if (*tiedpredIn==1){ 
		pairs[s] += 1/ww;
		conc[s] += 1/(2* ww);
	      }
	    }
	    else{
	      /*
		call pair concordant if p_i < p_j
	      */
	      pairs[s] += 1/ww;
	      if (pred[i + s * (*N)] < pred[j + s * (*N)]) {
		conc[s] += 1/ww;
	      }
	    }
	  }
	}
      }
    }
    AUC[s]=conc[s]/pairs[s];
    /* lasttime=times[s]; */
  }
}
