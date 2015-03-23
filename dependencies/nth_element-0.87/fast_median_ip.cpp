/* fast_median_ip.cpp
 * Ver 0.87
 * Peter H. Li 2013 FreeBSD License
 * See fast_median_ip.m for documentation.
 *
 * This uses undocumented MathWorks internals and should be used at your
 * own risk!
 */
#include "mex.h"
#include "unshare.hpp"
#include "fast_median_lib.cpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  // Check inputs
  if (nrhs != 1) {
    mexErrMsgIdAndTxt("Numerical:fast_median:nrhs", "Arguments should be the matrix of columns and the rank of the desired element");
  }
  if (!mxIsNumeric(prhs[0])) {
    mexErrMsgIdAndTxt("Numerical:fast_median:prhs", "Input argument must be a numeric matrix.");
  }

  // Unshare input array if necessary, then modify inplace.  UNDOCUMENTED MEX CALLS!!
  plhs[0] = (mxArray*) prhs[0];
  mxUnshareArray(plhs[0], true);
  plhs[0] = run_fast_median(plhs[0]);
}
