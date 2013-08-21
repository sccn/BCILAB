/*
 * CORE_resetintegrator.c
 * ----------------------
 * The syntax for this function is described in the M file of
 * the same name.
 *
 * This is a MEX-file for MATLAB.  It was compiled using:
 * 
 *      COREmex CORE_resetintegrator.c
 * 
 * samar@cs.stanford.edu, 02/20/06
 */

#include "matrix.h"
#include "mex.h"
#include <math.h>
#include "CORE_library.h"


/* ********************************************************** */
/* ********************* MEX GATEWAY  *********************** */
/* ********************************************************** */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  int N, mrows, ncols, j;
  double *input, *output;
  mxLogical *reset;

  /* ************* ACCESS MATLAB VARIABLES ************* */
  mrows = mxGetM(prhs[0]);   ncols = mxGetN(prhs[0]);
  N = mxGetNumberOfElements(prhs[0]);
  input = (double *) mxGetPr(prhs[0]);
  reset = mxGetLogicals(prhs[1]);

  /* Create an array to hold the output. */
  plhs[0] = mxCreateDoubleMatrix(mrows, ncols, mxREAL);
  output = mxGetPr(plhs[0]);

  /* ********************* COMPUTE ********************* */
  
  if (reset[0])  output[0] = input[0];  else  output[0] = 0;
  for (j = 1; j < N; j++) {
    if (reset[j])  output[j] = output[j-1] + input[j];
    else           output[j] = 0;
  }

  /* ********************* CLEAN UP ******************** */

}
