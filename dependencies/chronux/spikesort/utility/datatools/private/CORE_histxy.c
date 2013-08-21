/*
 * CORE_histxy.c
 * --------------------
 * The syntax for this function is described in the M file of
 * the same name.
 *
 * This is a MEX-file for MATLAB.  It was compiled using:
 * 
 *      COREmex CORE_histxy.c
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
  int outrows, outcols, count, i;
  double *rowdata, *coldata, *output;

  /* ************* ACCESS MATLAB VARIABLES ************* */
  count = mxGetM(prhs[0]);
  outrows = (int) mxGetScalar(prhs[3]);
  outcols = (int) mxGetScalar(prhs[2]);

  /* Get a handle to the input arrays. */
  rowdata = (double *) mxGetPr(prhs[1]);
  coldata = (double *) mxGetPr(prhs[0]);

  /* Create an array to hold the output. */
  plhs[0] = mxCreateNumericMatrix(outrows, outcols, mxDOUBLE_CLASS, mxREAL);
  output = (double *) mxGetPr(plhs[0]);

  /* ********************* COMPUTE ********************* */
  count_2D(rowdata, coldata, count, output, outrows, outcols);

  /* ********************* CLEAN UP ******************** */

}
