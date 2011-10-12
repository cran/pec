#include <R.h>
void cindex(double *C,
	    double *conc,
	    double *pairs,
	    int *tindex,
	    double *Y,
	    int *status,
	    double *times,
	    double *weight_i,
	    double *weight_j,
	    double *pred,
	    int *N,
	    int *NT,
	    int *tiedpredIn,
	    int *tiedoutcomeIn,
	    int *tiedmatchIn,
	    int *cens_model){
  int i,j,s;
  double wi, wj, lasttime=0;
  for (s=0; s<(*NT);s++) {
    conc[s]=0;
    pairs[s]=0;
    for (i=0;i<(*N);i++){
      /*
	for usuable pairs the smaller time must be uncensored
      */
      if (Y[i]<=times[s] && status[i]==1){ 
	for (j=i+1;j<*N;j++){
	  if (*cens_model==0){
	    /*
	      marginal censoring survival weights:
	      G(T_i-) for i
	      G(T_i)  for j
	    */
	    wi = weight_i[(tindex[i]-1)];
	    wj = weight_j[(tindex[i]-1)];
	  }
	  else{
	    /*
	      conditional censoring survival weights:
	      G(T_i-|X_i) for i
	      G(T_i|X_j)  for j
	    */
	    wi = weight_i[(tindex[i]-1)];
	    wj = weight_j[(j + (tindex[i]-1) * (*N))];
	  }
	  /*
	    pair unusuable if any weight==0
	  */
	  if (wj>0 && wi>0){
	    /*
	      rare case: same outcome and same prediction
	      count as concordant pair when
	      tiedmatchIn == TRUE
	    */
	    if (*tiedmatchIn==1 && (Y[i]==Y[j] && status[j]==1 && (pred[i + s * (*N)] == pred[j + s * (*N)]))){
	      pairs[s] += 1/(wi * wj);
	      conc[s] += 1/(wi * wj);
	    }
	    else{
	      /*
		if tiedoutcomeIn==0 call pairs with tied outcome unusuable,
		unless Y_j was censored, since then the uncensored Y_j will
		be greater than Y_i 
	      */
	      if (*tiedoutcomeIn==1 || (Y[i]!=Y[j] || status[j]==0)){
		if (pred[i + s * (*N)] == pred[j + s * (*N)]) {
		  /*
		    call pair unusuable if same predictions
		  */
		  if (*tiedpredIn==1){ 
		    pairs[s] += 1/(wi * wj);
		    conc[s] += 1/(2* wi * wj);
		  }
		}
		else{
		  /*
		    call pair concordant if p_i < p_j
		  */
		  pairs[s] += 1/(wi * wj);
		  if (pred[i + s * (*N)] < pred[j + s * (*N)]) {
		    conc[s] += 1/(wi * wj);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    C[s]=conc[s]/pairs[s];
    lasttime=times[s];
  }
}


