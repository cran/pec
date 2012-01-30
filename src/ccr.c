/*
  COMPETING RISK 
  Wolbers et al.
  conditional censoring survival weights:
  WAij=G(T_i-|X_i)G(T_i|X_j)
  WBij=G(T_i-|X_i)G(T_j-|X_j) 

  NOTE: there is one weight for each person in weight.i
  and one row for each person in weight.j
  we need the value at time T[i]
*/

#include <R.h>
void ccr(double *C,
	 double *concA,
	 double *pairsA,
	 double *concB,
	 double *pairsB,
	 int *tindex,
	 double *T,
	 int *Delta,
	 int *D,
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
  double Aij, Bij, WAij=1, WBij=1, weightedA, weightedB ,lasttime=0, weightedConcPairs,weightedPairs;
  for (s=0; s<(*NT);s++) {
    concA[s]=0;  /* count concordant pairs with (T[i]<T[j], D[i]=1,T[i]) */
    concB[s]=0;  /* count concordant pairs with (T[i]>=T[j], D[j]=2,T[i]) */
    pairsA[s]=0; /* count pairs with (T[i]<T[j], D[i]=1,T[i]) */
    pairsB[s]=0; /* count pairs with (T[i]>=T[j], D[j]=2,T[i]) */
    weightedConcPairs=0; /* weighted concordant pairs */
    weightedPairs=0;     /* weighted pairs */
    for (i=0;i<(*N);i++){
      /* for all pairs one of the times must be uncensored */
      if (T[i]<=times[s] && Delta[i]==1 && D[i]==1){
	/* Rprintf("\n\ni=%d\n",i+1); */
	for (j=0;j<*N;j++){
	  if (j!=i){
	    Aij=0;
	    Bij=0;
	    /* extract the weights */
	    if (*cens_model==0){
	      WAij = weight_i[i] * weight_j[(tindex[i]-1)];
	    }
	    else{
	      WAij = weight_i[i] * weight_j[(j + (tindex[i]-1) * (*N))];
	    }
	    WBij = weight_i[i] * weight_i[j];
	    /* time j must either be greater than time i */
	    if (T[i]<T[j] || (T[j]==T[i] && Delta[j]==0)){
	       /* Rprintf("case A i=%d\tj=%d\tT.i=%1.0f\tT.j=%1.0f\t\tD.i=%d\tD.j=%d\n",i+1,j+1,T[i],T[j],D[i],D[j]);  */
	      Aij=1;
	      weightedA=Aij/WAij;
	      Bij=0;
	      weightedB=0;
	    }
	    else{ /* or competing risk */
	      if(Delta[j]==1 && D[j]!=1){
		 /* Rprintf("case B i=%d\tj=%d\tT.i=%1.0f\tT.j=%1.0f\t\tD.i=%d\tD.j=%d\n",i+1,j+1,T[i],T[j],D[i],D[j]);  */
		Bij=1;
		weightedB=Bij/WBij;
		Aij=0;
		weightedA=0;
	      }
	      else{
		 /* Rprintf("case C i=%d\tj=%d\tT.i=%1.0f\tT.j=%1.0f\t\tD.i=%d\tD.j=%d\n",i+1,j+1,T[i],T[j],D[i],D[j]);  */
		Aij=0;
		Bij=0;
		weightedA=0;
		weightedB=0;
	      }
	    }
	    /* Rprintf("i=%d\tj=%d\tT[i]=%1.2f\tT[j]=%1.2f\tD[j]=%d\tpred[i]=%1.2f\tpred[j]=%1.2f\tAij=%1.2f\tBij=%1.2f\n",i,j,T[i],T[j],D[j],pred[i + s * (*N)],pred[j + s * (*N)],Aij,Bij);  */
	    /* Rprintf("s=%d\ti=%d\tj=%d\tAij=%1.2f\tBij=%1.2f\tpairsA.s=%1.2f\tpairsB.s=%1.2f\n",s,i,j,Aij,Bij,pairsA[s],pairsB[s]);  */
	    /* count pairs */
	    pairsA[s] += Aij;
	    pairsB[s] += Bij;
	    weightedPairs += (weightedA+weightedB);
	    /*

	      concordant pairs

	    */
	    if (pred[i + s * (*N)] > pred[j + s * (*N)]) {
	      concA[s] +=Aij;
	      concB[s] +=Bij;
	      weightedConcPairs += (weightedA+weightedB);
	      /* Rprintf("Concordant: pred.i=%1.2f\tpred.j=%1.2f\n",pred[i + s * (*N)],pred[j + s * (*N)]); */
	    } 
	    /*
	      pairs with equal predictions
	      count 1/2 or nothing
	    */
	    if (pred[i + s * (*N)] == pred[j + s * (*N)]) {
	      if (*tiedpredIn==1 ){ 
		concA[s] += (1/2) * Aij;
		concB[s] += (1/2) * Bij;
		weightedConcPairs += (1/2) * (weightedA+weightedB);
	      }
	    }
	  }
	}
      }
    }
    /* C[s]=(concA[s]+concB[s])/(pairsA[s]+pairsB[s]); */
    C[s]=weightedConcPairs/weightedPairs;
    lasttime=times[s];
  }
}
