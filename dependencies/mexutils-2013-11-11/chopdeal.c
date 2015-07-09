#include "mex.h"

#define SIZE_ARR 1
#define DATA_ARR 0

/* See chopdeal.m for documentation of this function. */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    int datalen,numsizes,k,s,p;
    double *indata_d, *outdata_d;
    float *indata_f, *outdata_f;
    char *indata_c, *outdata_c;
    double *sizes;
    int cursize;
    mxArray *tmp;
    if (nrhs != 2)
        mexErrMsgIdAndTxt("chopdeal:incorrect_input_argument_number","2 input arguments required."); 
    if (mxGetClassID(prhs[SIZE_ARR]) != mxDOUBLE_CLASS)
        mexErrMsgIdAndTxt("chopdeal:chopping_must_be_double","The chop-size array must be of type double."); 
    numsizes = mxGetNumberOfElements(prhs[SIZE_ARR]);
    datalen = mxGetNumberOfElements(prhs[DATA_ARR]);
    if (nlhs != numsizes)
        mexErrMsgIdAndTxt("chopdeal:incorrect_output_argument_number","The number of elements in the sizes array (%i) must equal the number of output arguments (%i).",numsizes,nlhs); 
    if (numsizes) {
        sizes = (double*)mxGetData(prhs[SIZE_ARR]);
        switch (mxGetClassID(prhs[DATA_ARR])) {
            case mxDOUBLE_CLASS:
                indata_d = (double*)mxGetData(prhs[DATA_ARR]);
                p = 0;
                for (s=0;s<numsizes;s++) {
                    cursize = sizes[s];
                    if (cursize<0)
                        mexErrMsgIdAndTxt("chopdeal:sizes_must_be_nonnegative","Sizes must be nonnegative."); 
                    plhs[s] = mxCreateDoubleMatrix(1,cursize,mxREAL);
                    if (cursize) {
                        outdata_d = (double*)mxGetData(plhs[s]);
                        for (k=0;k<cursize;k++,p++) {
                            if (p>=datalen)
                                mexErrMsgIdAndTxt("chopdeal:index_exceeds_matrix_dimensions","Index exceeds matrix dimensions."); 
                            *outdata_d++ = *indata_d++;
                        }
                    }
                }
                break;
            case mxSINGLE_CLASS:
                indata_f = (float*)mxGetData(prhs[DATA_ARR]);
                p = 0;
                for (s=0;s<numsizes;s++) {
                    cursize = sizes[s];
                    if (cursize<0)
                        mexErrMsgIdAndTxt("chopdeal:sizes_must_be_nonnegative","Sizes must be nonnegative."); 
                    plhs[s] = mxCreateNumericMatrix(1,cursize,mxSINGLE_CLASS,mxREAL);
                    if (cursize) {
                        outdata_f = (float*)mxGetData(plhs[s]);
                        for (k=0;k<cursize;k++,p++) {
                            if (p>=datalen)
                                mexErrMsgIdAndTxt("chopdeal:index_exceeds_matrix_dimensions","Index exceeds matrix dimensions."); 
                            *outdata_f++ = *indata_f++;
                        }
                    }
                }
                break;            
            case mxUINT8_CLASS:
                indata_c = (char*)mxGetData(prhs[DATA_ARR]);
                p = 0;
                for (s=0;s<numsizes;s++) {
                    cursize = sizes[s];
                    if (cursize<0)
                        mexErrMsgIdAndTxt("chopdeal:sizes_must_be_nonnegative","Sizes must be nonnegative."); 
                    plhs[s] = mxCreateNumericMatrix(1,cursize,mxUINT8_CLASS,mxREAL);
                    if (cursize) {
                        outdata_c = (char*)mxGetData(plhs[s]);
                        for (k=0;k<cursize;k++,p++) {
                            if (p>=datalen)
                                mexErrMsgIdAndTxt("chopdeal:index_exceeds_matrix_dimensions","Index exceeds matrix dimensions."); 
                            *outdata_c++ = *indata_c++;
                        }
                    }
                }
                break;            
            case mxCELL_CLASS:
                p = 0;
                for (s=0;s<numsizes;s++) {
                    cursize = sizes[s];
                    if (cursize<0)
                        mexErrMsgIdAndTxt("chopdeal:sizes_must_be_nonnegative","Sizes must be nonnegative."); 
                    plhs[s] = mxCreateCellMatrix(1,cursize);
                    if (cursize) {
                        for (k=0;k<cursize;k++,p++) {
                            if (p>=datalen)
                                mexErrMsgIdAndTxt("chopdeal:index_exceeds_matrix_dimensions","Index exceeds matrix dimensions.");
                            mxSetCell(plhs[s],k,mxDuplicateArray(mxGetCell(prhs[DATA_ARR],p)));
                        }
                    }
                }
                break;
            default:
                mexErrMsgIdAndTxt("chopdeal:unsuppored_class","Input array class not supported: %s",mxGetClassName(prhs[0]));
        }
    }
}
