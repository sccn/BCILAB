/*
 * fast_kmeans_step.c - fast pairwise distance matrix
 *
 * This MEX file works, but comments and documentation in particular are 
 *  a little slim ...
 * 
 *  Syntax:
 *    [assigns,bestdists,clustersizes,centroids] = 
 *                   fast_kmeans(waves,centroids,normsqr_waves,chunk)
 * 
 *    ##### NOTE ##### Waves and centroids must be of size:
 *            waves:     d rows, m cols 
 *            centroids: d rows, n cols
 * 
 *  (needs a LAPACK library (libmwlapack.lib (Borland) , 
 *  msvc_libmwlapack.lib (MSVC), lcc_libmwlapack.lib (LCC)) to 
 *   compile, e.g.:   mex fast_kmeans_step.c borland_libmwlapack.lib     )
 * 
 * This is a MEX-file for MATLAB.
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"
#include "matrix.h"

#define MIN(a,b)  (((a) < (b)) ? (a) : (b))



/* Prototypes */
double *normsqr_vects(double *vects, int mcols, int drows);
void minus2xty(double* output, double *vects1, double *vects2,
               int mcols, int ncols, int drows);
void psst_init(int mrows, int ncols);
void psst_cleanup(void);


/* Variables */
static double *psst_chunk = NULL;
static int psst_m = 0;
static int psst_n = 0;


/* Function: normsqr(double *, int, int)
 * ------------------------------------------------------------------
 * Compute norms (along the columns) of a list of vectors.
 */
double *normsqr_vects(double *vects, int mcols, int drows)
{
  int m, d;
  double *normvec, temp;

  normvec = (double *) mxMalloc (mcols * sizeof(double));

  for (m = 0; m < mcols; m++) {
    normvec[m] = 0;
    for (d = 0; d < drows; d++) {
      temp = vects[m*drows + d];
      normvec[m] += temp*temp;
    }
  }
  return normvec;
}

/* Function: minus2xty(double *,double *,double *,PTR_dgemm,int,int,int)
 * ------------------------------------------------------------------
 * Subtracts 2*vects1(:,i)*vects2(:,i)' from each element of output.
 */
void minus2xty(double* output, double *vects1, double *vects2,
               int mcols, int ncols, int drows)
{
  int i, numel;
  double minus2 = -2.0, one = 1.0;
  char chn = 'N', cht = 'T';

  dgemm(&cht, &chn, &mcols, &ncols, &drows, &minus2,
	  vects1, &drows, vects2, &drows, &one, output, &mcols);

  numel = ncols*mcols;
  for (i = 0; i < numel; i++) {    // Fix roundoff error ...
    if (output[i] < 0) {output[i] = 0;}
  }
}


/* Function: psst_init()
 * ------------------------------------------------------------------
 * Allocates persistent memory if necessary.
 */
void psst_init(int mrows, int ncols)
{
  if (psst_chunk != NULL) {
    if ((psst_m != mrows) || (psst_n != ncols))
      { mxFree(psst_chunk);  psst_chunk = NULL; psst_m = 0; psst_n = 0;}
  }
  if (psst_chunk == NULL) {
    psst_chunk = (double *) mxMalloc(mrows*ncols*sizeof(double));
    psst_m = mrows;  psst_n = ncols;
    mexMakeMemoryPersistent(psst_chunk);
    mexAtExit(psst_cleanup);
  }
}



/* Function: psst_cleanup()
 * ------------------------------------------------------------------
 * Cleans up persistent memory..
 */
void psst_cleanup(void) {
  mxFree(psst_chunk);
  psst_chunk = NULL;    psst_m = 0;   psst_n = 0;
}


/* ********************************************************** */
/* ********************* MEX GATEWAY  *********************** */
/* ********************************************************** */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  int mcols, ncols, drows, K, m, n, d, singleton, idx;
  int chunksize, numchunks, chunk, offset, chunklen;
  double *data1, *data2, *normsqr1, *normsqr2, value;
  double *distsq, *bestdists, *clustsizes, *centroids;
  unsigned long *assigns;
  int junk;  

  mcols = mxGetN(prhs[0]);  drows = mxGetM(prhs[0]);  ncols = mxGetN(prhs[1]);
  data1     = (double *) mxGetPr(prhs[0]);
  data2     = (double *) mxGetPr(prhs[1]);
  normsqr1  = (double *) mxGetPr(prhs[2]);
  chunksize = (int) mxGetScalar(prhs[3]);
  chunksize = MIN(chunksize, mcols);

  /* Create an array to hold the outputs (if needed). */
  psst_init(chunksize, ncols);
  plhs[0] = mxCreateNumericMatrix(mcols, 1, mxUINT32_CLASS, mxREAL);
  plhs[1] = mxCreateNumericMatrix(mcols, 1, mxDOUBLE_CLASS, mxREAL);
  plhs[2] = mxCreateNumericMatrix(ncols, 1, mxDOUBLE_CLASS, mxREAL);
  plhs[3] = mxCreateNumericMatrix(drows, ncols, mxDOUBLE_CLASS, mxREAL);
  assigns    = (unsigned long *) mxGetPr(plhs[0]);
  bestdists  = (double *) mxGetPr(plhs[1]);
  clustsizes = (double *) mxGetPr(plhs[2]);
  centroids  = (double *) mxGetPr(plhs[3]);

  distsq   = psst_chunk;  
  //plhs[4] = mxCreateNumericMatrix(chunksize,ncols, mxDOUBLE_CLASS, mxREAL);
  //distsq = (double *) mxGetPr(plhs[4]);

  /* ********************* COMPUTE ********************* */
  normsqr2 = normsqr_vects(data2, ncols, drows);  

  // Distances are computed as follows:
  //     dist(x,y)^2 = (x-y)'(x-y) = x'x + y'y - 2 x'y
  // but we process waveforms in chunks (of chunksize each)
  // at once to save memory.

  numchunks = (int)ceil((double)mcols/(double)chunksize);
  for (chunk = 0; chunk < numchunks; chunk++) {
    offset = chunk * chunksize;
    chunklen = MIN(chunksize, mcols - offset);

    // First, fill in the x'x and y'y terms ...
    for (n = 0; n < ncols; n++) {
      for (m = 0; m < chunklen; m++) {
	distsq[m + n*chunklen] = normsqr1[m + offset] + normsqr2[n];
      }
    }
    
    // Next comes the -2x'y part ...
    minus2xty(distsq, &data1[offset*drows], data2, chunklen, ncols, drows);
    
    // Make cluster assignments
    for (m = 0; m < chunklen; m++) {  idx = m+offset;
      bestdists[idx] = distsq[m];  assigns[idx] = 0;
      for (n = 0; n < ncols; n++) {
	if (distsq[m+n*chunklen] < bestdists[idx])
	  {  bestdists[idx] = distsq[m+n*chunklen];  assigns[idx] = n; }
      }
      clustsizes[assigns[idx]]++;
      assigns[idx] += 1;  // convert to 1-indexed for Matlab's sake
    }
  }
  
  // Recompute centroids
  for (m = 0; m < mcols; m++) {
    for (d = 0; d < drows; d++)
      { centroids[(assigns[m]-1)*drows + d] += data1[m*drows + d]; }
  }

  for (n = 0; n < ncols; n++) {
    for (d = 0; d < drows; d++)
      { centroids[n*drows + d] /= clustsizes[n]; }
  }

  /* ********************* CLEAN UP ******************** */
  mxFree(normsqr2);
}
