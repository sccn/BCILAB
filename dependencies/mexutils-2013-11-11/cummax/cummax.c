#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  double *dp, *dend, dmax;
  float *fp, *fend, fmax;
  if (nrhs != 1)
      mexErrMsgIdAndTxt("cummax:one_input_required","This function has one input argument.");
  if (nlhs > 1)
      mexErrMsgIdAndTxt("cummax:one_output_required","This function has one output argument.");
  plhs[0] = mxDuplicateArray(prhs[0]);
  if (mxIsEmpty(plhs[0]))
    return;
  switch (mxGetClassID(plhs[0])) {
      case mxDOUBLE_CLASS:
          dp = mxGetData(plhs[0]);
          dend = dp + mxGetNumberOfElements(plhs[0]);
          dmax = *dp++;
          while (dp != dend)
              if (*dp < dmax)
                  *dp++ = dmax;
              else
                  dmax = *dp++;
          break;
      case mxSINGLE_CLASS:
          fp = mxGetData(plhs[0]);
          fend = fp + mxGetNumberOfElements(plhs[0]);
          fmax = *fp++;
          while (fp != fend)
              if (*fp < fmax)
                  *fp++ = fmax;
              else
                  fmax = *fp++;
          break;
      default:
          mexErrMsgIdAndTxt("cummax:unsuppored_class","Input array class not supported: %s",mxGetClassName(plhs[0]));
  }  
}
