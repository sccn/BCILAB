#include "ssa_matlab.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*******************************************************\
     MATLAB
     ycom = signal_int(alpha, mcom, Z, Zdyn, t0, postproc)
    \*******************************************************/
    
    /*** Input variables ***/
    int p, m, n, t0;
    int Zd;
    double *Z, *alpha, *dmcom;
    int ncom, *mcom, postproc;
    
    /*** Output variables ***/
    int dims[3];
    double *ycom;
    
    /*** Temporary variables ***/
    int i;

    /*** Retrieve data ***/
    m           = mxGetM(prhs[0]);
    n           = mxGetN(prhs[0]);
    alpha       = mxGetPr(prhs[0]);

    /*** Retrieve model ***/
    ncom        = mxGetNumberOfElements(prhs[1]) - 1;
    dmcom       = mxGetPr(prhs[1]);
    mcom        = (int *)mxCalloc(ncom, sizeof(int));
    for(i=0; i<ncom; i++) mcom[i] = (int)(dmcom[i+1] - dmcom[i]);
    get_M(prhs + 2, &p, &Z, &Zd);
    t0          = (int)mxGetScalar(prhs[4]);
    
    /*** Retrieve analysis settings ***/
    postproc    = (int)mxGetScalar(prhs[5]);
    
    /*** Allocate output variables ***/
    dims[0] = p; dims[1] = n; dims[2] = ncom;
    ycom    = mxGetPr(plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
    
    /*** Extract signal ***/
    signal(p, m, n, ncom, mcom, Zd ? Z + p*p*(t0-1) : Z, Zd, alpha, ycom);
    
    /*** Cleanup ***/
    mxFree(mcom);
    cleanup();

    /*** Postprocessing ***/
    if(postproc && p == 1) {
        permute_data(plhs[0], 1, n, ncom);
        dims[0] = ncom; dims[1] = n;
        mxSetDimensions(plhs[0], dims, 2);
    }
}
