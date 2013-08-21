/**
 *  char[] osc_address_new('ipaddr / hostname', port_num)
 *
 */

#include "mex.h"
#include "lo/lo.h"

char maddr[32];
char port[32];
char host[256];

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  int err;
  int hslen;
  lo_address d;

  if (nrhs != 2) {
    mexErrMsgTxt("Please specify host name/address and port");
    return;
  }
  if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments.");
    return;
  }

  if(! mxIsChar(prhs[0])) {
    mexErrMsgTxt("Expecting a character array in the first argument.");
    return;
  }

  if(! mxIsNumeric(prhs[1])) {
    mexErrMsgTxt("Expecting a numeric scalar in the second argument.");
    return;
  }

  hslen = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
  if(hslen > 255) {
    mexErrMsgTxt("Maximum host name string length is 255 characters.");
    return;
  }
  err = mxGetString(prhs[0], host, hslen);

  if(err != 0) {
    mexErrMsgTxt("Error reading host name string.");
    return;
  }

  sprintf(port, "%d", (int)(mxGetScalar(prhs[1])));
  d = lo_address_new(host, port);

  if(! d) {
    mexErrMsgTxt(lo_address_errstr(d));
    return;
  }

  sprintf(maddr, "osc_address:%x", (void*) d);
  plhs[0] = mxCreateString(maddr);

}
