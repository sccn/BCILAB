/*
 * CORE_leadingedges.c
 * --------------------
 * The syntax for this function is described in the M file of
 * the same name.
 *
 * This is a MEX-file for MATLAB.  It was compiled using:
 * 
 *      COREmex CORE_leadingedges.c
 * 
 * samar@cs.stanford.edu, 01/20/06
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
  int mrows, ncols, i, j, jcol;
  double *inputD64;
  unsigned char *inputU8;
  mxLogical *output;

  /* ************* ACCESS MATLAB VARIABLES ************* */
  mrows = mxGetM(prhs[0]);  ncols = mxGetN(prhs[0]);

  /* Create an array to hold the output. */
  plhs[0] = mxCreateLogicalMatrix(mrows, ncols);
  output = mxGetLogicals(plhs[0]);

  /* ********************* COMPUTE ********************* */

  switch (mxGetClassID(prhs[0])) {
    case (mxDOUBLE_CLASS):
      inputD64 = (double *) mxGetPr(prhs[0]);
      for (j = 0; j < ncols; j++) {
	jcol = j*mrows;
	output[0+jcol] = (inputD64[0+jcol] != 0);
	for (i = 1; i < mrows; i++) {
	  output[i+jcol] = inputD64[i+jcol] && !inputD64[(i-1)+jcol];
	}
      }
      break;
    case(mxUINT8_CLASS):
    case(mxLOGICAL_CLASS):
      inputU8 = (unsigned char *) mxGetData(prhs[0]);
      for (j = 0; j < ncols; j++) {
	jcol = j*mrows;
	output[0+jcol] = (inputU8[0+jcol] != 0);
	for (i = 1; i < mrows; i++) {
	  output[i+jcol] = inputU8[i+jcol] && !inputU8[(i-1)+jcol];
	}
      }
      break;
    default:
      mexErrMsgTxt("Invalid Data Type.\n");
  }
  

  /* ********************* CLEAN UP ******************** */

}
