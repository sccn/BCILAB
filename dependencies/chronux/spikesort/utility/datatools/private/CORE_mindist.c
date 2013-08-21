/*
 * CORE_mindist.c
 * --------------------
 * The syntax for this function is described in the M file of
 * the same name.
 *
 * This is a MEX-file for MATLAB.  It was compiled using:
 *      COREmex CORE_mindist.c
 * 
 * samar@cs.stanford.edu, 09/12/05
 */

#include "mex.h"
#include "matrix.h"
#include "CORE_library.h"


/* ********************************************************** */
/* ********************* MEX GATEWAY  *********************** */
/* ********************************************************** */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  int mrows, nrows, dcols;
  double *vectsX, *normsqrX, *vectsY, *normsqrY, *yinds, *dists;

  /* ************* ACCESS MATLAB VARIABLES ************* */
  /* Get all of the flags first -- no error checking, so they better be ok */

  /* Next get pointer to vector lists */
  vectsX  = (double *) mxGetPr(prhs[0]);   vectsY = (double *)mxGetPr(prhs[1]);
  mrows = mxGetM(prhs[0]);   dcols = mxGetN(prhs[0]);    nrows = mxGetM(prhs[1]);

  /* Compute vector norms */
  normsqrX = (double *) mxMalloc(mrows * sizeof(double));
  norm_squared_rows(vectsX, mrows, dcols, normsqrX);

  normsqrY = (double *) mxMalloc(nrows * sizeof(double));
  norm_squared_rows(vectsY, nrows, dcols, normsqrY);


  /* ************** ALLOCATE OUTPUT MATRIX ************* */
  plhs[0] = mxCreateNumericMatrix(mrows, 1, mxDOUBLE_CLASS, mxREAL);
  yinds = (double *) mxGetPr(plhs[0]);
  plhs[1] = mxCreateNumericMatrix(mrows, 1, mxDOUBLE_CLASS, mxREAL);
  dists = (double *) mxGetPr(plhs[1]);

  /* ********************* COMPUTE ********************* */
  minimum_distances(vectsX, normsqrX, vectsY, normsqrY, 
                      mrows, nrows, dcols, yinds, dists);

  /* ********************* CLEAN UP ******************** */
  mxFree(normsqrX); 
  mxFree(normsqrY);
}
