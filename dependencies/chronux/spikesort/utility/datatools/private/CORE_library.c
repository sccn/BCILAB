/*
 * FILE: CORE_library.c
 * --------------------
 * Common computational routines for Matlab CORE_function_ MEX files.  
 * Function descriptions can be found in CORE_library.h
 * 
 * samar@cs.stanford.edu, 01/18/06
 */

#include <math.h>
#include <limits.h>
#include <float.h>
#include "CORE_library.h"



/* ------------------------------------------------------------------------
 * Function: count_2D
 * ------------------------------------------------------------------------ */
void count_2D(double *rowinds, double *colinds, int N, 
              double *counts, int maxrow, int maxcol)
{
  int i;

  for (i = 0; i < N; i++) {
    counts[(int)(rowinds[i]-1) + maxrow*((int)colinds[i]-1)]++;
  }
}
/* ------------------------------------------------------------------------ */


/* ------------------------------------------------------------------------
 * Function: count_diffs
 * ------------------------------------------------------------------------ */

void count_diffs(double *x, int Nx, double *y, int Ny, double *diffs, int maxdiff)
{
  int curr_y, lo_x, hi_x, curr_x, offset;

  lo_x = 0;  hi_x = 0;
  for (curr_y = 0; curr_y < Ny; curr_y++) {
    while ((y[curr_y]-x[lo_x]) > maxdiff) {   // move lo index till < maxdiff before y[curr]
      lo_x++;
      if (lo_x > Nx-1) {curr_y = Ny; break;}
    }
    if (hi_x < lo_x) hi_x = lo_x;             // catch up hi index if it was passed
    if (hi_x <= Nx-1) {                     
      while ((x[hi_x]-y[curr_y]) <= maxdiff) {// move hi index until > maxdiff past y[curr]
	hi_x++;
	if (hi_x > Nx-1)  break;
      }
    }

    offset = -(int)y[curr_y] + maxdiff;       // <- no need to do this in the loop below   
    for (curr_x=lo_x; curr_x<hi_x; curr_x++) {// now count each value btw the indices
      diffs[(int)x[curr_x]+offset]++;
    }
  }
}
/* ------------------------------------------------------------------------ */


/* ------------------------------------------------------------------------
 * Function: count_in_cols
 * ------------------------------------------------------------------------ */
void count_in_cols(double *matrix, int M, int N, double *counts, int maxval)
{
  int i, j;

  for (i = 0; i < N; i++) {     // for each col,
    for (j = 0; j < M; j++) {     // go through the rows, ...
      counts[((int)matrix[j + i*M]-1) + i*maxval]++;  // and tabulate
    }
  }
}
/* ------------------------------------------------------------------------ */


/* ------------------------------------------------------------------------
 * Function: find_minmax
 * ------------------------------------------------------------------------ */
void find_minmax(double *input, int M, 
                 double *minval, int *minind, double *maxval, int *maxind)
{
  int i;

  for (i = 0; i < M; i++) {
    if (!IsNaN(input[i])) {   // init min/max to first non-NAN element
      *minval =  input[i];    *minind = i;   
      *maxval =  input[i];    *maxind = i;
      break;
    }
  }
  if (i == M) {  // special case: if all NaNs
    *minval = GetNaN();  *minind = 0;   *maxval = GetNaN();  *maxind = 0;
  }

  for (; i < M; i++) {   // go through keeping track of min/max
    if      ((input[i] < *minval) && (!IsNaN(input[i]))) 
        {*minval = input[i];  *minind = i;}
    else if ((input[i] > *maxval) && (!IsNaN(input[i]))) 
        {*maxval = input[i];  *maxind = i;}
  }
}
/* ------------------------------------------------------------------------ */


/* ------------------------------------------------------------------------
 * Function: k_centers_approx
 * ------------------------------------------------------------------------ */
void k_centers_approx(double *vects, int M, int N, int K, double *centers)
{

}
// ------------------------------------------------------------------------


/* ------------------------------------------------------------------------
 * Function: k_means
 * ------------------------------------------------------------------------ */
void k_means(double *vects, int M, int N, int K, double *centroids)
{

}
/* ------------------------------------------------------------------------ */


/* ------------------------------------------------------------------------
 * Function: minimum_column_indices
 * ------------------------------------------------------------------------ */
void minimum_column_indices(double *matrix, int M, int N, double *mininds)
{
  int i, j, soFarInd;
  double soFar, nextTry;

  for (i = 0; i < M; i++) {
    soFar = matrix[i];  soFarInd = 0;
    for (j = 0; j < N; j++) {
      nextTry = matrix[i + j*M];
      if (nextTry < soFar) {soFar = nextTry;  soFarInd = j;}
    }
    mininds[i] = soFarInd + 1;  // change to one index
  }
}
/* ------------------------------------------------------------------------ */


/* ------------------------------------------------------------------------
 * Function: minimum_distances
 * ------------------------------------------------------------------------ */
void minimum_distances(double *vectsX, double *normsqrX, double *vectsY, double *normsqrY,
                           int M, int N, int D, double *Yinds, double *distsqr)
{
  int i, j, k, soFarInd;
  double soFar, nextTry, temp;
  
  // (Distances are computed as in pairdist_squared_rows below)
  for (i = 0; i < M; i++) {
    soFar = DBL_MAX;  soFarInd = 0;
    for (j = 0; j < N; j++) {
      temp = 0.0;
      for (k = 0; k < D; k++) {
	temp += (vectsX[i+k*M] * vectsY[j+k*N]);
      }
      nextTry = normsqrX[i] + normsqrY[j] - 2*temp;
      if (nextTry < soFar) {soFar = nextTry;  soFarInd = j;}
    }
    distsqr[i] = soFar;  Yinds[i] = soFarInd + 1;
  }

  // alternative distance computation: (x-y)'(x-y)
  /*nextTry = 0.0;
  for (k = 0; k < D; k++) {
    temp = (vectsX[i+k*M] - vectsY[j+k*N]);
    nextTry += temp * temp;
  }*/
}
/* ------------------------------------------------------------------------ */


/* ------------------------------------------------------------------------
 * Function: norm_squared_rows
 * ------------------------------------------------------------------------ */
void norm_squared_rows(double *matrix, int M, int N, double *normsqr)
{
  int i, j;
  double temp;

  for (i = 0; i < M; i++) {
    normsqr[i] = 0;
    for (j = 0; j < N; j++) {
      temp = matrix[i + j*M];   // avoid matrix access 2x in next line
      normsqr[i] += temp*temp;
    }
  }
}
// ------------------------------------------------------------------------


/* ------------------------------------------------------------------------
 * Function: pairdist_squared_rows
 * ------------------------------------------------------------------------ */
void pairdist_squared_rows(double *vectsX, double *normsqrX, double *vectsY, double *normsqrY,
                           double *dists, int M, int N, int D)
{
  int m,n;
  char chn = 'N', cht = 'T';
  double minus2 = -2.0, one = 1.0;

  /* Distances are computed as follows:
   *    dist(x,y)^2 = (x-y)'(x-y) = x'x + y'y - 2 x'y
   * The 2nd formulation requires fewer flops than the 1st. */
  for (n = 0; n < N; n++) {
    for (m = 0; m < M; m++) {
      dists[m + n*M] = normsqrX[m] + normsqrY[n];  // write x'x + y'y
    }
  }  
  dgemm(&chn, &cht, &M, &N, &D, &minus2, vectsX, &M, vectsY, &N, &one, dists, &M);  // then add -2x'y

  if ((vectsX == vectsY) && (M == N)) {  // Self-distance; force diagonal to 0
    for (n = 0; n < N; n++) {  dists[n + n*N] = 0;  }
  }
}
/* ------------------------------------------------------------------------ */


/* ------------------------------------------------------------------------
 * Function: pairdist_squared_rows_safe
 * ------------------------------------------------------------------------ */
void pairdist_squared_rows_safe(double *vectsX, double *vectsY, double *dists, 
                                int M, int N, int D)
{
  int m, n, d;
  double temp;

  /* Distances are computed as follows:
   *    dist(x,y)^2 = (x-y)'(x-y)
   * giving dist(x,x)^2 == 0 w/o the roundoff error of the prev fcn. */
  for (n = 0; n < N; n++) {
    for (m = 0; m < N; m++) {
      dists[m + n*M] = 0;
      for (d = 0; d < D; d++) {
	temp = vectsX[m + d*M] - vectsY[n + d*N];
	dists[m + n*M] = temp*temp;
      }
    }
  }
}
/* ------------------------------------------------------------------------ */


/* ------------------------------------------------------------------------
 * Function: sqrt_matrix_inplace
 * ------------------------------------------------------------------------ */
void sqrt_matrix_inplace(double *matrix, int M, int N)
{
  int i,numel;
  
  numel = N*M;
  for (i = 0; i < numel; i++) {
    matrix[i] = sqrt(matrix[i]);
  }
}
/* ------------------------------------------------------------------------ */


