/*
 * CORE_minmax.c
 * --------------------
 * The syntax for this function is described in the M file of
 * the same name.
 *
 * This is a MEX-file for MATLAB.  It was compiled using:
 * 
 *      COREmex CORE_minmax.c
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
  int N, minind, maxind;
  double *input, *xmin, *xmax, *minindout, *maxindout;

  /* ************* ACCESS MATLAB VARIABLES ************* */
  N = mxGetNumberOfElements(prhs[0]);

  /* Get a handle to the input array. */
  input = (double *) mxGetPr(prhs[0]);

  /* Create arrays to hold the outputs. */
  plhs[0] = mxCreateDoubleScalar(0);  // xmin
  xmin = mxGetPr(plhs[0]);
  plhs[1] = mxCreateDoubleScalar(0);  // xmax
  xmax = mxGetPr(plhs[1]);
  plhs[2] = mxCreateDoubleScalar(0);  // minind
  minindout = mxGetPr(plhs[2]);
  plhs[3] = mxCreateDoubleScalar(0);  // maxind
  maxindout = mxGetPr(plhs[3]);

  /* ********************* COMPUTE ********************* */
  find_minmax(input, N, xmin, &minind, xmax, &maxind);

  /* ******************** CLEAN UP ********************* */
  *minindout = (double) minind + 1; // 0-index int => 1-index dbl
  *maxindout = (double) maxind + 1; //    (for Matlab's sake)
}


