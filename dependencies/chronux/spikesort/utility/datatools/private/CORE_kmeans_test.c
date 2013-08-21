/*
 * fast_kmeans_step.c - fast pairwise distance matrix
 *
 * Warning: this MEX file works, but is still being updated.
 *    Comments and documentation in particular are a little
 *    slim ...
 * 
 *  Syntax:
 *    [assigns,bestdists,clustersizes,centroids] = 
 *                   fast_kmeans(waves,centroids,normsqr_waves,chunk)
 * 
 *    ##### NOTE ##### Waves and centroids must be of size:
 *            waves:     d rows, m cols 
 *            centroids: d rows, n cols
 * 
 *  (needs libmwlapack.lib (Borland) or msvc_libmwlapack.lib (MSVC) to
 *   compile:   mex fast_kmeans_step.c libmwlapack.lib     )
 * 
 * This is a MEX-file for MATLAB.
 */

#include <float.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "CORE_library.h"


/* ********************************************************** */
/* ********************* MEX GATEWAY  *********************** */
/* ********************************************************** */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  int mrows, nrows, dcols, K, m, n, d;
  double *data1, *data2, *normsqr1, *normsqr2;
  double *assigns, *bestdists, *clustsizes, *centroids;

  mrows = mxGetM(prhs[0]);  nrows = mxGetM(prhs[1]);  dcols = mxGetN(prhs[1]);
  data1     = (double *) mxGetPr(prhs[0]);
  data2     = (double *) mxGetPr(prhs[1]);
  normsqr1  = (double *) mxGetPr(prhs[2]);

  /* Create an array to hold the outputs (if needed). */
  plhs[0] = mxCreateNumericMatrix(mrows, 1, mxUINT32_CLASS, mxREAL);
  plhs[1] = mxCreateNumericMatrix(mrows, 1, mxDOUBLE_CLASS, mxREAL);
  plhs[2] = mxCreateNumericMatrix(nrows, 1, mxDOUBLE_CLASS, mxREAL);
  plhs[3] = mxCreateNumericMatrix(nrows, dcols, mxDOUBLE_CLASS, mxREAL);
  assigns    = (double *) mxGetPr(plhs[0]);
  bestdists  = (double *) mxGetPr(plhs[1]);
  clustsizes = (double *) mxGetPr(plhs[2]);
  centroids  = (double *) mxGetPr(plhs[3]);

  /* ********************* COMPUTE ********************* */
  norm_squared_rows(data2, nrows, dcols, normsqr2);  
  minimum_distances(data1, normsqr1, data2, normsqr2, 
                      mrows, nrows, dcols, assigns, bestdists);

  // Count membership
  for (m = 0; m < mrows; m++) {  clustsizes[(int)assigns[m]]++;  }

  // Recompute centroids
  for (m = 0; m < mrows; m++) {
    for (d = 0; d < dcols; d++)
      { centroids[((int)assigns[m]-1) + d*mrows] += data1[m + d*mrows]; }
  }

  for (n = 0; n < nrows; n++) {
    for (d = 0; d < dcols; d++)
      { centroids[n + d*nrows] /= clustsizes[n]; }
  }

  /* ********************* CLEAN UP ******************** */
  mxFree(normsqr2);
}
