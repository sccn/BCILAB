/*
 * fast_interface_energy.c - Matlab mex file
 *
 * This is a MEX-file for MATLAB.
 */

#include "mex.h"
#include "matrix.h"

/* ********************************************************** */
/* ********************* MEX GATEWAY  *********************** */
/* ********************************************************** */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  int mrows, ncols, i, numel;
  double *input, *output, energy, scale;


  /* ************ MATLAB ARGUMENT CHECKING ************* */
  if (nrhs != 2)  {mexErrMsgTxt("Function takes two inputs."); }
  if (nlhs > 1) { mexErrMsgTxt("Function only returns one output."); }
  if (mxGetNumberOfDimensions(prhs[0]) != 2)
    { mexErrMsgTxt("First input must be a 2-D matrix."); }
  if ((mxGetM(prhs[1]) != 1) || (mxGetN(prhs[1]) != 1))
  { mexErrMsgTxt("Second input must be a scalar."); }

  mrows = mxGetM(prhs[0]);  ncols = mxGetN(prhs[0]);

  /* Get a handle to the inputs. */
  input = (double *) mxGetPr(prhs[0]);
  scale = -1/mxGetScalar(prhs[1]);

  /* Create an scalar to hold the output. */
  plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
  output = (double *) mxGetPr(plhs[0]);

  /* ********************* COMPUTE ********************* */
  energy = 0;
  numel = mrows * ncols;
  for (i = 0; i < numel; i++) {
    energy += exp(input[i] * scale);
  }
  *output = energy;

  /* ********************* CLEAN UP ******************** */
}
