#include "ssa_matlab.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /****************************************************************************************************************\
     MATLAB
     [snlogL rrN uuD] = loglik_grad_int(y, H, Hd, Z, Zd, T, Td, R, Rd, Q, Qd, c, cd, a1, P1, tol0, tolP, inv_method)
    \****************************************************************************************************************/
    
    /*** Input variables ***/
    int p, m, r, n;
    double *y;
    int *miss, *pmiss, *imiss;
    int Hd, Zd, Td, Rd, Qd, cd, sta;
    double *H, *Z, *T, *R, *Q, *c, *a1, *P1, *P1_inf;
    double tol0, tolP;
    int inv_method;
    
    /*** Output variables ***/
    int dims[3];
    double *snlogL, *rrN, *uuD;
    
    /*** Buffer variables ***/
    double *v, *invF, *K, *L;
    int d, *Fns;
    
    /*** Retrieve data ***/
    get_data(prhs[0], &p, &n, (int *)0, &y, &miss);
    
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

    /*** Retrieve analysis settings ***/
    tol0        = mxGetScalar(prhs[15]);
    tolP        = mxGetScalar(prhs[16]);
    inv_method  = (int)mxGetScalar(prhs[17]);

    /*** Truncate Z and H w.r.t. missing data ***/
    if(miss) truncate_yHZ(p, m, n, 1, &y, miss, &pmiss, &imiss, &H, &Hd, &Z, &Zd);
    else { pmiss = (int *)0; imiss = (int *)0; }
    
    /*** Allocate buffer ***/
    Fns     = (int *)mxCalloc(n, sizeof(int));
    v       = (double *)mxCalloc(p*n, sizeof(double));
    invF    = (double *)mxCalloc(p*p*n, sizeof(double));
    K       = (double *)mxCalloc(m*p*n, sizeof(double));
    L       = (double *)mxCalloc(m*m*n, sizeof(double));
    
    /*** Allocate output variables ***/
    snlogL  = mxGetPr(plhs[0] = mxCreateDoubleScalar(0));
    dims[0] = m; dims[1] = m; dims[2] = n;
    rrN     = mxGetPr(plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
    dims[0] = p; dims[1] = p; dims[2] = n;
    uuD     = mxGetPr(plhs[2] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
    
    /*** Kalman filter ***/
    kalman(p, m, r, n, 1, y, (int *)0, false, pmiss, sta, H, Hd, Z, Zd, T, Td, R, Rd, Q, Qd, c, cd, a1, P1, P1_inf,
        (double *)0, (double *)0, (double *)0, &d, Fns, v, (double *)0, invF, (double *)0, K, (double *)0,
        L, (double *)0, snlogL, (double *)0, (double *)0, (double *)0, (double *)0, tol0, tolP, inv_method);
    
    /*** Smoothing cumulant ***/
    smocum(p, m, n, (int *)0, false, pmiss, Z, Zd, T, Td, d, Fns, v, invF, K, L, uuD, rrN);
    
    /*** Cleanup ***/
    mxFree(Fns);
    mxFree(v);
    mxFree(invF);
    mxFree(K);
    mxFree(L);
    cleanup();
}
