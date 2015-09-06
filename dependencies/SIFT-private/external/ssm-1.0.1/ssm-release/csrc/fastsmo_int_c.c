#include "ssa_matlab.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /********************************************************************************************************************\
     MATLAB
     [alphahat epshat etahat] = fastsmo_int(y, H, Hd, Z, Zd, T, Td, R, Rd, Q, Qd, c, cd, a1, P1,
                                            tol0, tolP, inv_method, postproc)
     [epshat etahat] = fastsmo_int(...)
     alphahat = fastsmo_int(...)
    \********************************************************************************************************************/
    
    /*** Input variables ***/
    int p, m, r, n, N;
    double *y;
    int *miss, *pmiss, *imiss;
    int Hd, Zd, Td, Rd, Qd, cd, sta;
    double *H, *Z, *T, *R, *Q, *c, *a1, *P1, *P1_inf;
    double tol0, tolP;
    int inv_method, postproc;
    
    /*** Output variables ***/
    int dims[3];
    double *alphahat, *epshat, *etahat;
    
    /*** Buffer variables ***/
    double *v, *invF, *F2, *K, *L, *L1, *QRt, *RQRt;
    int d, *Fns;
    
    /*** Retrieve data ***/
    get_data(prhs[0], &p, &n, &N, &y, &miss);
    
    /*** Retrieve model ***/
    get_M(prhs + 1, &p, &H, &Hd);
    get_M(prhs + 3, (int *)0, &Z, &Zd);
    get_M(prhs + 5, &m, &T, &Td);
    get_M(prhs + 7, (int *)0, &R, &Rd);
    get_M(prhs + 9, &r, &Q, &Qd);
    get_M(prhs + 11, (int *)0, &c, &cd);
    get_M(prhs + 13, (int *)0, &a1, (int *)0);
    get_P1(prhs[14], (int *)0, &P1, &P1_inf);
    sta     = !Hd && !Zd && !Td && !Rd && !Qd;
    
    /*** Retrieve analysis setting ***/
    tol0        = mxGetScalar(prhs[15]);
    tolP        = mxGetScalar(prhs[16]);
    inv_method  = (int)mxGetScalar(prhs[17]);
    postproc    = (int)mxGetScalar(prhs[18]);

    /*** Truncate Z and H w.r.t. missing data ***/
    if(miss) truncate_yHZ(p, m, n, N, &y, miss, &pmiss, &imiss, &H, &Hd, &Z, &Zd);
    else { pmiss = (int *)0; imiss = (int *)0; }
    
    /*** Allocate output variables ***/
    switch(nlhs) {
        case 1:
            dims[0]     = m; dims[1] = N; dims[2] = n;
            alphahat    = mxGetPr(plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
            epshat      = (double *)0;
            etahat      = (double *)0;
            break;
        case 2:
            alphahat    = (double *)0;
            dims[0]     = p; dims[1] = N; dims[2] = n;
            epshat      = mxGetPr(plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
            dims[0]     = r; dims[1] = N; dims[2] = n;
            etahat      = mxGetPr(plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
            break;
        default:
            dims[0]     = m; dims[1] = N; dims[2] = n;
            alphahat    = mxGetPr(plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
            dims[0]     = p; dims[1] = N; dims[2] = n;
            epshat      = mxGetPr(plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
            dims[0]     = r; dims[1] = N; dims[2] = n;
            etahat      = mxGetPr(plhs[2] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
            break;
    }
    
    /*** Allocate workspace ***/
    Fns     = (int *)mxCalloc(n, sizeof(int));
    v       = (double *)mxCalloc(p*N*n, sizeof(double));
    invF    = (double *)mxCalloc(p*p*n, sizeof(double));
    F2      = (double *)mxCalloc(p*p*n, sizeof(double));
    K       = (double *)mxCalloc(m*p*n, sizeof(double));
    L       = (double *)mxCalloc(m*m*n, sizeof(double));
    L1      = (double *)mxCalloc(m*m*n, sizeof(double));
    QRt     = (double *)((Rd || Qd) ? mxCalloc(r*m*n, sizeof(double)) : mxCalloc(r*m, sizeof(double)));
    RQRt    = (double *)((Rd || Qd) ? mxCalloc(m*m*n, sizeof(double)) : mxCalloc(m*m, sizeof(double)));

    /*** Kalman filter ***/
    kalman(p, m, r, n, N, y, (int *)0, false, pmiss, sta, H, Hd, Z, Zd, T, Td, R, Rd, Q, Qd, c, cd, a1, P1, P1_inf,
        (double *)0, (double *)0, (double *)0, &d, Fns, v, (double *)0, invF, F2, K, (double *)0,
        L, L1, (double *)0, (double *)0, (double *)0, QRt, RQRt, tol0, tolP, inv_method);
    
    /*** Fast smoother ***/
    fastsmo(p, m, r, n, N, (int *)0, false, pmiss, H, Hd, Z, Zd, T, Td, R, Rd, c, cd, a1, P1, P1_inf,
        d, Fns, v, invF, F2, K, L, L1, QRt, RQRt, Rd || Qd, (double *)0, (double *)0,
        alphahat, epshat, etahat);
    
    /*** Cleanup ***/
    mxFree(Fns);
    mxFree(v);
    mxFree(invF);
    mxFree(F2);
    mxFree(K);
    mxFree(L);
    mxFree(L1);
    mxFree(QRt);
    mxFree(RQRt);
    cleanup();

    /*** Postprocessing ***/
    if(postproc) {
        switch(nlhs) {
            case 1:
                permute_data(plhs[0], m, N, n);
                break;
            case 2:
                permute_data(plhs[0], p, N, n);
                permute_data(plhs[1], r, N, n);
                break;
            default:
                permute_data(plhs[0], m, N, n);
                permute_data(plhs[1], p, N, n);
                permute_data(plhs[2], r, N, n);
                break;
        }
    }
}
