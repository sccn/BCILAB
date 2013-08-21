/**
 *  void osc_server_free('osc_server:........')
 *
 */

#include "mex.h"
#include "lo/lo.h"

char maddr[32];

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  int err;
  lo_server s;
  int mlen;
  if (nrhs != 1) {
    mexErrMsgTxt("Expecting one argument");
    return;
  }
  if (nlhs > 0) {
    mexErrMsgTxt("Too many output arguments.");
    return;
  }

  if(! mxIsChar(prhs[0])) {
    mexErrMsgTxt("Expecting a character array in the first argument.");
    return;
  }

  mlen = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
  if(mlen > 32) {
    mexErrMsgTxt("Input string is too long.");
    return;
  }
  err = mxGetString(prhs[0], maddr, mlen);

  if(err != 0) {
    mexErrMsgTxt("Error reading input string.");
    return;
  }

  err = sscanf(maddr, "osc_server:%x", &s);

  if(err < 0 || (! s)) {
    mexErrMsgTxt("Error scanning input string.");
    return;
  }

  lo_server_free(s);
  
}
