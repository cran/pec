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
  double wi, wj, ww, lasttime=0;
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
	    wi = weight_i[i];
	    wj = weight_j[(tindex[i]-1)];
	  }
	  else{
	    /*
	      conditional censoring survival weights:
	      G(T_i-|X_i) for i
	      G(T_i|X_j)  for j

	      NOTE: there is one weight for each person in weight.i
	            and one row for each person in weight.j
		    we need the value at time Y[i]
	    */
	    /* wi = weight_i[(tindex[i]-1)]; */
	    wi = weight_i[i];
	    wj = weight_j[(j + (tindex[i]-1) * (*N))];
	  }
	  ww = (wi * wj);
	  /* Rprintf("i=%d\twi=%1.8f\n",i,wi); */
	  /* if ((1/wi)>1000) Rprintf("Large wi=%1.8f\n",wi); */
	  /* if ((1/ww)>100) Rprintf("Yi=%1.2f\tYj=%1.2f\tt=%1.2f\tLarge wj=%1.8f\twi=%1.8f\n",Y[i],Y[j],times[s],wi,wj);  */
	  /* if ((1/ww)>((*N)*(*N))) ww=1/((*N)*(*N)); */
	  /* if ((1/ww)>trunc[1]) {	   */
	  /* ww=1/trunc[1]; */
	  /* } */
	  /* if ((1/ww)>(*N)) { */
	  /* Rprintf("Warning: truncated weights"); */
	  /* Rprintf("Yi=%1.2f\tYj=%1.2f\twj=%1.8f\twi=%1.8f\n",Y[i],Y[j],wj,wi); */
	  /* ww=1/((double)(*N)); */
	  /* } */
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
	      pairs[s] += (1/ww);
	      conc[s] += (1/ww);
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
      }
    }
    C[s]=conc[s]/pairs[s];
    lasttime=times[s];
  }
}


