#include "bsb_common.h"
        
/*
 * bsb_close(hDLL, hConn, hUSB_WRITE, hCLOSE_DRIVER_ASYNC, pBuffer);
 * Backend function for BioSemi interface; do not use directly.
 * 
 * In:
 *   hDLL   : handle to the driver DLL to be closed
 *   hConn  : handle to the driver connection to be closed
 *   hUSB_WRITE            : handle to a driver function
 *   hCLOSE_DRIVER_ASYNC   : handle to a driver function
 *   pBuffer : pointer to the circular buffer
 *   
 */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] ) 
{
    uintptr_t *pTmp;  /* temp pointer */
    void *hDLL;     /* handle to the DLL */
    void *hConn;    /* handle to the driver connection */
    void *pBuffer;  /* pointer to the ring buffer */
    /* function handles */
    USB_WRITE_t USB_WRITE;
    CLOSE_DRIVER_ASYNC_t CLOSE_DRIVER_ASYNC;
    char zeros[64] = {0};
    
    if (nrhs != 5)
        mexErrMsgTxt("Five input argument required."); 
    if (1 != (mxGetNumberOfElements(prhs[0]) * mxGetNumberOfElements(prhs[1]) * mxGetNumberOfElements(prhs[2]) * mxGetNumberOfElements(prhs[3]) * mxGetNumberOfElements(prhs[4])))
        mexErrMsgTxt("All inputs must be scalars.");
    
    /* get DLL handle */
    pTmp = (uintptr_t*)mxGetData(prhs[0]);
    hDLL = (void*)*pTmp;

    /* get driver connection handle */
    pTmp = (uintptr_t*)mxGetData(prhs[1]);
    hConn = (void*)*pTmp;
    
    /* get handshake disable function */
    pTmp = (uintptr_t*)mxGetData(prhs[2]);
    USB_WRITE = (USB_WRITE_t)*pTmp;

    /* get driver close function */
    pTmp = (uintptr_t*)mxGetData(prhs[3]);
    CLOSE_DRIVER_ASYNC = (CLOSE_DRIVER_ASYNC_t)*pTmp;

    /* get buffer pointer */
    pTmp = (uintptr_t*)mxGetData(prhs[4]);
    pBuffer = (void*)*pTmp;

    /* close driver connection */  
    if (!USB_WRITE(hConn,zeros))
        mexWarnMsgIdAndTxt("BioSemi:HandshakeProblem","The handshake while shutting down the BioSemi driver gave an error.");
    if (!CLOSE_DRIVER_ASYNC(hConn))
        mexWarnMsgIdAndTxt("BioSemi:CloseProblem","Closing the BioSemi driver gave an error.");   
    
    /* close the DLL */
    FREE_LIBRARY(hDLL);

    /* free the ring buffer */
    free(pBuffer);
}

