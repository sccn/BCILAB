/*
 * CORE_histxt.c
 * --------------------
 * The syntax for this function is described in the M file of
 * the same name.
 *
 * This is a MEX-file for MATLAB.  It was compiled using:
 * 
 *      COREmex CORE_histxt.c
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
  int inrows, outrows, ncols, i, j;
  double *data, *output;

  /* ************* ACCESS MATLAB VARIABLES ************* */
  ncols = mxGetN(prhs[0]);
  inrows = mxGetM(prhs[0]);
  outrows = (int) mxGetScalar(prhs[1]);

  /* Get a handle to the input arrays. */
  data = (double *) mxGetPr(prhs[0]);

  /* Create an array to hold the output. */
  plhs[0] = mxCreateNumericMatrix(outrows, ncols, mxDOUBLE_CLASS, mxREAL);
  output = (double *) mxGetPr(plhs[0]);


  /* ********************* COMPUTE ********************* */
  count_in_cols(data, inrows, ncols, output, outrows);

  /* ********************* CLEAN UP ******************** */

}
