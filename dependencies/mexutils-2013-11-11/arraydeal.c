#include "mex.h"
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    int numel,k;
    void *data;
    if (nrhs != 1)
        mexErrMsgIdAndTxt("arraydeal:incorrect_input_argument_number","1 input argument required."); 
    numel = mxGetNumberOfElements(prhs[0]);
    if (nlhs != numel)
        mexErrMsgIdAndTxt("arraydeal:incorrect_output_argument_number","The number of elements in the input array must equal the number of output arguments."); 
    if (numel) {
        data = (void*)mxGetData(prhs[0]);
        switch (mxGetClassID(prhs[0])) {
            case mxDOUBLE_CLASS:
                for (k=0;k<numel;k++)
                    plhs[k] = mxCreateDoubleScalar(((double*)data)[k]);
                break;
            case mxSINGLE_CLASS:
                for (k=0;k<numel;k++) {
                    plhs[k] = mxCreateNumericMatrix(1,1,mxSINGLE_CLASS,mxREAL);
                    *(float*)mxGetData(plhs[k]) = ((float*)data)[k];
                }
                break;            
            default:
                mexErrMsgIdAndTxt("arraydeal:unsuppored_class","Input array class not supported: %s",mxGetClassName(prhs[0]));
        }
    }
}
