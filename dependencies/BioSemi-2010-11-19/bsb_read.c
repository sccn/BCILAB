#include "bsb_common.h"

/*
 * last_ptr = bsb_read(hDLL, hConn, hREAD_POINTER, pBuffer, nbchan, last_ptr);
 * Backend function for BioSemi interface; do not use directly.
 * 
 * In:
 *   hDLL       : handle to the driver DLL to be closed
 *   hConn      : handle to the driver connection to be closed
 *   hREAD_POINTER : handle to a driver function
 *   pBuffer    : pointer to the circular buffer
 *   nbchan     : number of channels
 *   last_ptr   : last byte offset into the buffer
 *
 * Out:
 *   RawData   : raw data block (signed integers, to be scaled by 1/8192 uV)
 *   last_ptr  : updated last byte offset into the buffer (needed for the next call)
 *
 */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] ) 
{
    uintptr_t *pTmp;  /* temp pointer */
    int32_T *pTmpI;   /* temp pointer2 */
    void *hDLL;       /* handle to the DLL */
    void *hConn;      /* handle to the driver connection */
    char *pBuffer;    /* pointer to the ring buffer */
    READ_POINTER_t READ_POINTER;
    int32_T last_ptr;  /* last buffer offset */
    int32_T cur_ptr;   /* pointer past the current sample */
    int32_T nbchan;    /* number of channels */    
    int32_T samplesize; /* size of a complete n-channel sample, in bytes */
    int32_T buffersize; /* size of the entire buffer, in bytes */
    
    mwSize msize[] = {0,0};  /* matrix size */
    mwSize scalar[] = {1,1}; /* matrix size */
    
    if (nrhs != 6)
        mexErrMsgTxt("Six input arguments required."); 
    if (1 != (mxGetNumberOfElements(prhs[0]) * mxGetNumberOfElements(prhs[1]) * mxGetNumberOfElements(prhs[2]) * mxGetNumberOfElements(prhs[3]) * mxGetNumberOfElements(prhs[4]) * mxGetNumberOfElements(prhs[5])))
        mexErrMsgTxt("All inputs must be scalars.");
    if (nlhs != 2)
        mexErrMsgTxt("Two output arguments required."); 
    
    /* get driver connection handle */
    pTmp = (uintptr_t*)mxGetData(prhs[1]);
    hConn = (void*)*pTmp;    
    /* get read pointer function */
    pTmp = (uintptr_t*)mxGetData(prhs[2]);
    READ_POINTER = (READ_POINTER_t)*pTmp;
    /* get buffer pointer */
    pTmp = (uintptr_t*)mxGetData(prhs[3]);
    pBuffer = (char*)*pTmp;
    /* get channel number */
    pTmp = (int32_T*)mxGetData(prhs[4]);
    nbchan = *pTmp;
    /* get last buffer offset */
    pTmp = (int32_T*)mxGetData(prhs[5]);
    last_ptr = *pTmp;
    
    /* compute a few sizes */
    samplesize = nbchan*4;
    buffersize = samplesize * BUFFER_SAMPLES;
        
    /* get current buffer offset */    
    if (!READ_POINTER(hConn, &cur_ptr)) {
        mexWarnMsgIdAndTxt("BioSemi:ReadProblem","Reading the updated buffer pointer gave an error.");
        cur_ptr = last_ptr;
    }
    
    /* forget about incomplete sample data */
    cur_ptr = cur_ptr  - cur_ptr % samplesize;
    if (cur_ptr < 0)
        cur_ptr = cur_ptr + buffersize;
        
    if (cur_ptr == last_ptr) {
        /* no new sample came in: create empty int32 data block */
        plhs[0] = mxCreateNumericArray(2,msize,mxINT32_CLASS,mxREAL);
    } else {
        msize[0] = nbchan;
        if (cur_ptr > last_ptr) {        
            /* sequential read: return intermediate part */
            msize[1] = (cur_ptr - last_ptr) / samplesize;
            plhs[0] = mxCreateNumericArray(2, msize, mxINT32_CLASS, mxREAL);
            memcpy((char*)mxGetData(plhs[0]), &pBuffer[last_ptr], cur_ptr-last_ptr);
        } else {
            /* wrap-around read: concatenate two parts */
            msize[1] = (cur_ptr + buffersize - last_ptr) / samplesize;
            plhs[0] = mxCreateNumericArray(2, msize, mxINT32_CLASS, mxREAL);
            memcpy((char*)mxGetData(plhs[0]), &pBuffer[last_ptr], buffersize-last_ptr);
            memcpy((char*)mxGetData(plhs[0])+buffersize-last_ptr, pBuffer, msize[0]*msize[1]*4-(buffersize-last_ptr));
        }
    }
    /* assign last_ptr output */
    plhs[1] = mxCreateNumericArray(2,scalar,mxINT32_CLASS,mxREAL); *(int32_T*)mxGetData(plhs[1]) = (int32_T)cur_ptr; 
}
