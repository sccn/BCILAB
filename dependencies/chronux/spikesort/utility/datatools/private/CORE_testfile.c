/*
 * CORE_testfile.c
 * --------------------
 * The syntax for this function is described in the M file of
 * the same name.
 *
 * This is a MEX-file for MATLAB.  It was compiled using:
 * 
 *      mex CORE_testfile.c
 * 
 * samar@cs.stanford.edu, 07/12/05
 */

#include "mex.h"
#include "matrix.h"
#include <math.h>

//#define DEBUG

/* ********************************************************** */
/* ********************* MEX GATEWAY  *********************** */
/* ********************************************************** */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  int mrows, ncols, i, j, temp, temp2;
  double *input1, *input2;
  mxLogical *output;

  /* ************* ACCESS MATLAB VARIABLES ************* */
  //  mrows = mxGetM(prhs[0]);  ncols = mxGetN(prhs[0]);

  /* Get a handle to the input array. */
  input1 = (double *) mxGetPr(prhs[0]);
  input2 = (double *) mxGetPr(prhs[1]);

  /* Create an array to hold the output. */
  plhs[0] = mxCreateLogicalScalar(true);
  output = (double *) mxGetLogicals(plhs[0]);


  /* ********************* COMPUTE ********************* */
  *output = false;
  if (*input1 > *input2) {  *output = true;  }

  /* ********************* CLEAN UP ******************** */

}
