/* nth_element_ip.cpp
 * Ver 0.86
 * Peter H. Li 2012 FreeBSD License 
 * See nth_element_ip.m for documentation.
 *
 * This uses undocumented MathWorks internals and should be used at your
 * own risk!
 */
#include "mex.h"
#include "unshare.hpp"
#include "nth_element_lib.cpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  // Check inputs
  if (nrhs != 2) {
    mexErrMsgIdAndTxt("Numerical:nth_element:nrhs", "Arguments should be the matrix of columns and the rank of the desired element");
  }
  if (!mxIsNumeric(prhs[0])) {
    mexErrMsgIdAndTxt("nth_element:prhs", "First argument must be a numeric matrix.");
  }
  if (!mxIsNumeric(prhs[1]) || mxGetNumberOfDimensions(prhs[1]) != 2 || mxGetM(prhs[1]) != 1 || mxGetN(prhs[1]) != 1) {
    mexErrMsgIdAndTxt("nth_element:prhs", "Second argument must be a scalar.");
  }
  const mwSize nrows = mxGetM(prhs[0]);
  const mwSize ncols = mxGetN(prhs[0]);


  // Validate rank argument
  mwIndex rank = (mwIndex) mxGetScalar(prhs[1]);
  if (rank < 1) {
    mexErrMsgIdAndTxt("nth_element:prhs", "Rank cannot be less than 1.");
  }
  if (rank > nrows) {
    mexErrMsgIdAndTxt("nth_element:prhs", "Rank cannot be greater than the number of rows.");
  }

  // Convert matlab-style index (starts at 1) to C++ (starts at 0).
  rank--;
  
  
  // If user wants to also get the rearranged indices, have to do extra work
  unsigned int *indices;
  if (nlhs > 1) {
      plhs[1] = mxCreateNumericMatrix(nrows, ncols, mxUINT32_CLASS, mxREAL);
      indices = (unsigned int *) mxGetData(plhs[1]);
      for (mwIndex i = 0; i < ncols; ++i) {
          for (mwIndex j = 0; j < nrows; ++j) indices[i*nrows + j] = j;
      }
  }
  
  
  // Unshare input array if necessary, then modify inplace.  UNDOCUMENTED CALLS!!
  plhs[0] = (mxArray*) prhs[0];
  mxUnshareArray(plhs[0], true);
  if (nlhs > 1) {
      run_nth_element(plhs[0], indices, rank, ncols, nrows);
      for (mwIndex i = 0; i < ncols*nrows; ++i) ++indices[i]; // Fix for janky Matlab 1-indexing
  } else {
    run_nth_element(plhs[0], rank, ncols, nrows);
  }
}
