void brier_noinf(double *bs,
		 double *Y,
		 double *pred,
		 int *N)
{
  int i, j;
  double p, y, brier;
  
  for (j=0; j<*N; j++){
    p = pred[j];    /* prediction */
    for (i=0; i<*N; i++){
      y = Y[i];	/* observation */
      brier=(y-p)*(y-p);
      *bs += brier / (double) ((*N) * (*N));
    }
  }
}

