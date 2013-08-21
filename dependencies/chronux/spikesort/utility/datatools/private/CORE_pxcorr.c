/*
 * CORE_pxcorr.c
 * --------------------
 * The syntax for this function is described in the M file of
 * the same name.
 *
 * This is a MEX-file for MATLAB.  It was compiled using:
 * 
 *      COREmex CORE_pxcorr.c
 * 
 * samar@cs.stanford.edu, 07/31/05
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
  int lenx, leny, maxlag, a, b, c, bb, offset;
  double *x, *y, *output;

  /* ************* ACCESS MATLAB VARIABLES ************* */
  lenx = mxGetN(prhs[0]);   leny = mxGetN(prhs[1]);

  /* Get a handle to the input array. */
  x = (double *) mxGetPr(prhs[0]);
  y = (double *) mxGetPr(prhs[1]);
  maxlag = (int) mxGetScalar(prhs[2]);

  /* Create an array to hold the output. */
  plhs[0] = mxCreateNumericMatrix(1, 2*maxlag+1, mxDOUBLE_CLASS, mxREAL);
  output = (double *) mxGetPr(plhs[0]);


  /* ********************* COMPUTE ********************* */  
  count_diffs(x, lenx, y, leny, output, maxlag);

  /* ********************* CLEAN UP ******************** */

}
