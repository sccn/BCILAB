/* fast_median.cpp
 * Ver 0.87
 * Peter H. Li 2013 FreeBSD License
 * See fast_median.m for documentation.
 */
#include "mex.h"
#include "fast_median_lib.cpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  // Check inputs
  if (nrhs != 1) {
    mexErrMsgIdAndTxt("fast_median:nrhs", "Arguments should be the matrix of columns and the rank of the desired element");
  }
  if (!mxIsNumeric(prhs[0])) {
    mexErrMsgIdAndTxt("fast_median:prhs", "Input argument must be a numeric matrix.");
  }

  // Copy input array, pass to generic method
  mxArray *incopy = mxDuplicateArray(prhs[0]);
  plhs[0] = run_fast_median(incopy);
}
