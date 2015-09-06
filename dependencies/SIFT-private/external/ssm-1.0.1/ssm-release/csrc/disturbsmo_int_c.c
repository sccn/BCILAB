#include "ssa_matlab.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /***************************************************************************************************\
     MATLAB
     [epshat etahat epsvarhat etavarhat r N]
        = disturbsmo_int(y, H, Hd, Z, Zd, T, Td, R, Rd, Q, Qd, c, cd, a1, P1,
                            tol0, tolP, inv_method, postproc)
    \***************************************************************************************************/
    
    /*** Input variables ***/
    int p, m, rr, n, NN;
    double *y;
    int *miss, *pmiss, *imiss;
    int Hd, Zd, Td, Rd, Qd, cd, sta;
    double *H, *Z, *T, *R, *Q, *c, *a1, *P1, *P1_inf;
    double tol0, tolP;
    int inv_method, postproc;
    
    /*** Output variables ***/
    int dims[3];
    double *epshat, *etahat, *epsvarhat, *etavarhat, *r, *N;
    
    /*** Buffer variables ***/
    double *v, *invF, *K, *L, *RQ, *QRt;
    int d, *Fns;
    
    /*** Retrieve data ***/
    get_data(prhs[0], &p, &n, &NN, &y, &miss);
    
    /*** Retrieve model ***/
    get_M(prhs + 1, &p, &H, &Hd);
    get_M(prhs + 3, (int *)0, &Z, &Zd);
    get_M(prhs + 5, &m, &T, &Td);
    get_M(prhs + 7, (int *)0, &R, &Rd);
    get_M(prhs + 9, &rr, &Q, &Qd);
    get_M(prhs + 11, (int *)0, &c, &cd);
    get_M(prhs + 13, (int *)0, &a1, (int *)0);
    get_P1(prhs[14], (int *)0, &P1, &P1_inf);
    sta     = !Hd && !Zd && !Td && !Rd && !Qd;
    
    /*** Retrieve analysis settings ***/
    tol0        = mxGetScalar(prhs[15]);
    tolP        = mxGetScalar(prhs[16]);
    inv_method  = (int)mxGetScalar(prhs[17]);
    postproc    = (int)mxGetScalar(prhs[18]);
    
    /*** Truncate Z and H w.r.t. missing data ***/
    if(miss) truncate_yHZ(p, m, n, NN, &y, miss, &pmiss, &imiss, &H, &Hd, &Z, &Zd);
    else { pmiss = (int *)0; imiss = (int *)0; }
    
    /*** Allocate output variables ***/
    dims[0] = p; dims[1] = NN; dims[2] = n;
    epshat  = mxGetPr(plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
    dims[0] = rr; dims[1] = NN; dims[2] = n;
    etahat  = mxGetPr(plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
    dims[0] = p; dims[1] = p; dims[2] = n;
    epsvarhat   = mxGetPr(plhs[2] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
    dims[0] = rr; dims[1] = rr; dims[2] = n;
    etavarhat   = mxGetPr(plhs[3] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
    dims[0] = m; dims[1] = NN; dims[2] = n+1;
    r       = mxGetPr(plhs[4] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
    dims[0] = m; dims[1] = m; dims[2] = n+1;
    N       = mxGetPr(plhs[5] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
    
    /*** Allocate workspace ***/
    Fns     = (int *)mxCalloc(n, sizeof(int));
    v       = (double *)mxCalloc(p*NN*n, sizeof(double));
    invF    = (double *)mxCalloc(p*p*n, sizeof(double));
    K       = (double *)mxCalloc(m*p*n, sizeof(double));
    L       = (double *)mxCalloc(m*m*n, sizeof(double));
    RQ      = (double *)((Rd || Qd) ? mxCalloc(m*rr*n, sizeof(double)) : mxCalloc(m*rr, sizeof(double)));
    QRt     = (double *)((Rd || Qd) ? mxCalloc(rr*m*n, sizeof(double)) : mxCalloc(rr*m, sizeof(double)));

    /*** Kalman filter ***/
    kalman(p, m, rr, n, NN, y, (int *)0, false, pmiss, sta, H, Hd, Z, Zd, T, Td, R, Rd, Q, Qd, c, cd, a1, P1, P1_inf,
        (double *)0, (double *)0, (double *)0, &d, Fns, v, (double *)0, invF, (double *)0, K, (double *)0,
        L, (double *)0, (double *)0, (double *)0, RQ, QRt, (double *)0, tol0, tolP, inv_method);
    
    /*** Fast smoother ***/
    disturbsmo(p, m, rr, n, NN, (int *)0, false, pmiss, H, Hd, Z, Zd, T, Td, Q, Qd, d, Fns, v, invF,
        K, L, RQ, QRt, Rd || Qd, r, N, epshat, etahat, epsvarhat, etavarhat);
    
    /*** Cleanup ***/
    mxFree(Fns);
    mxFree(v);
    mxFree(invF);
    mxFree(K);
    mxFree(L);
    mxFree(RQ);
    mxFree(QRt);
    cleanup();
    
    /*** Postprocessing ***/
    if(postproc) {
        permute_data(plhs[0], p, NN, n);
        permute_data(plhs[1], rr, NN, n);
        permute_data(plhs[4], m, NN, n+1);
        if(p == 1) permute_data(plhs[2], 1, 1, n);
        if(rr == 1) permute_data(plhs[3], 1, 1, n);
        if(m == 1) permute_data(plhs[5], 1, 1, n+1);
    }
}
