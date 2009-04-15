void cindex(double *C,
	    int *tindex,
	    double *Y,
	    int *status,
	    double *weight,
	    double *weight_lag,
	    double *Z,
	    int *N,
	    double *tau,
	    int *cens_model){
  int i,j;
  double concord=0, comparison=0, wi, wj;
  for (i=0;i<(*N);i++){
    if (Y[i]<=*tau && status[i]==1){
      for (j=i+1;j<*N;j++){
	if (*cens_model==0){
	  wi = weight_lag[(tindex[i]-1)];
	  wj = weight[(tindex[i]-1)];
	  }
	else{
	  wi = weight_lag[(i + (tindex[i]-1) * (*N))];
	  wj = weight[(j + (tindex[i]-1) * (*N))];
	}
	if (wj>0 && wi>0){ 
	  comparison += 1/(wi * wj);
	  if (Z[i]<Z[j]){
	    concord += 1/(wi * wj);
	  } 
	}
      }
    }
  }
  C[0]=concord/comparison;
}



void hindex(double *C,
	    double *Y,
	    int *status,
	    double *Z,
	    int *N,
	    int *tau){
  int i,j;
  double concord=0, comparison=0;
  for (i=0;i<(*N-*tau);i++){
    for (j=i+1;j<*N;j++){
      if (status[i]==1){
	comparison++;
	if (Z[i]<Z[j])
	  concord++;
      }
    }
  }
  C[0]=concord/comparison;
}


/* void cindex(double *C, */ 
/* 	    int *status, */
/* 	    double *surv, */
/* 	    double *bisurv, */
/* 	    double *G, */
/* 	    double *GZ, */
/* 	    int *N){ */
  
/*   int i, j, s, x; */
/*   double w=0,v=0; */
  
/*   for (i=0;i<(*N);i++){ */
/*     if (status[i]==1){ */
/*       if (surv[i]>0) */
/* 	w += surv[i] / G[i]; */
/*       if (bisurv[i]>0) */
/* 	v += bisurv[i] / GZ[i]; */
/*     } */
/*     printf("G=%1.2f\tS=%1.2f\tSZ=%1.2f\tGZ=%1.2f\tv=%1.2f\tw=%1.2f\n",surv[i],G[i],bisurv[i],GZ[i],v,w);   */ 
/*   } */
/*   C[0]=v/w; */
/* } */
