#include "ssa_matlab.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*******************************************************************************\
     MATLAB
     [wt_a wt_alpha] = weights_int(n, t0, H, Hd, Z, Zd, T, Td, R, Rd, Q, Qd, P1,
                                    com, mcom, tol0, tolP, inv_method, postproc)
    \*******************************************************************************/
    
    /*** Input variables ***/
    int p, m, r, n, t0;
    int Hd, Zd, Td, Rd, Qd, sta;
    double *H, *Z, *T, *R, *Q, *P1, *P1_inf;
    double tol0, tolP;
    int com, inv_method, postproc;
    
    /*** Output variables ***/
    double *wt_a, *wt_alpha, *wt_fil, *wt_smo;
    
    /*** Retrieve model ***/
    n       = (int)mxGetScalar(prhs[0]);
    t0      = (int)mxGetScalar(prhs[1]);
    get_M(prhs + 2, &p, &H, &Hd);
    get_M(prhs + 4, (int *)0, &Z, &Zd);
    get_M(prhs + 6, &m, &T, &Td);
    get_M(prhs + 8, (int *)0, &R, &Rd);
    get_M(prhs + 10, &r, &Q, &Qd);
    get_P1(prhs[12], (int *)0, &P1, &P1_inf);
    sta     = !Hd && !Zd && !Td && !Rd && !Qd;
    
    /*** Retrieve analysis settings ***/
    com         = (int)mxGetScalar(prhs[13]);
    tol0        = mxGetScalar(prhs[15]);
    tolP        = mxGetScalar(prhs[16]);
    inv_method  = (int)mxGetScalar(prhs[17]);
    postproc    = (int)mxGetScalar(prhs[18]);

    if(com) {
        /*** Retrieve model component dimensions ***/
        int i, dims[3];
        int ncom        = mxGetNumberOfElements(prhs[14]) - 1;
        double *dmcom   = mxGetPr(prhs[14]);
        int *mcom       = (int *)mxCalloc(ncom, sizeof(int));
        for(i=0; i<ncom; i++) mcom[i] = (int)(dmcom[i+1] - dmcom[i]);
        
        /*** Allocate output variables ***/
        wt_a    = (double *)mxCalloc(m*p*n, sizeof(double));
        dims[0] = p; dims[1] = p*n; dims[2] = ncom;
        wt_fil  = mxGetPr(plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
        if(nlhs > 1) {
            wt_alpha    = (double *)mxCalloc(m*p*n, sizeof(double));
            wt_smo      = mxGetPr(plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
        }
        else wt_alpha   = (double *)0;
        
        /*** Weights ***/
        weights(p, m, r, n, t0, (int *)0, false, sta, H, Hd, Z, Zd, T, Td, R, Rd, Q, Qd, P1, P1_inf, wt_a, wt_alpha, tol0, tolP, inv_method);
        signal(p, m, p*n, ncom, mcom, Zd ? Z + p*p*(t0-1) : Z, false, wt_a, wt_fil);
        if(nlhs > 1) {
            signal(p, m, p*n, ncom, mcom, Zd ? Z + p*p*(t0-1) : Z, false, wt_alpha, wt_smo);
            mxFree(wt_alpha);
        }
        
        /*** Postprocessing ***/
        if(postproc && p == 1) {
            dims[0] = ncom; dims[1] = n;
            permute_data(plhs[0], 1, n, ncom);
            mxSetDimensions(plhs[0], dims, 2);
            if(nlhs > 1) {
                permute_data(plhs[1], 1, n, ncom);
                mxSetDimensions(plhs[1], dims, 2);
            }
        }

        /*** Cleanup ***/
        mxFree(mcom);
        mxFree(wt_a);
    }
    else {
        /*** Allocate output variables ***/
        wt_a        = mxGetPr(plhs[0] = mxCreateDoubleMatrix(m, p*n, mxREAL));
        wt_alpha    = (nlhs > 1) ? mxGetPr(plhs[1] = mxCreateDoubleMatrix(m, p*n, mxREAL)) : (double *)0;
        
        /*** Weights ***/
        weights(p, m, r, n, t0, (int *)0, false, sta, H, Hd, Z, Zd, T, Td, R, Rd, Q, Qd, P1, P1_inf, wt_a, wt_alpha, tol0, tolP, inv_method);
    }

    /*** Cleanup ***/
    cleanup();
}
