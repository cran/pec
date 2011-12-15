#include <math.h>
#include <R.h>
void SNull(double *time,
	   double *jumptimes,
	   double *elp,
	   double *S,
	   int *N,
	   int *NJ){
  int s,i;
  for (s=0; s<*NJ; s++){
    for (i=0; i<*N; i++){
      if (time[i]>=jumptimes[s])
	S[s]+=elp[i];
    }
  }
}
