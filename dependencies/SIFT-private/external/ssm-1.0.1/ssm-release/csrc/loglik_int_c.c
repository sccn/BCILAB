#include "ssa_matlab.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /**************************************************************************************************************\
     MATLAB
     [snlogL fvar] = loglik_int(y, H, Hd, Z, Zd, T, Td, R, Rd, Q, Qd, c, cd, a1, P1, tol0, tolP, inv_method)
    \**************************************************************************************************************/
    
    /*** Input variables ***/
    int p, m, r, n, N;
    double *y;
    int *miss, *pmiss, *imiss;
    int Hd, Zd, Td, Rd, Qd, cd, sta;
    double *H, *Z, *T, *R, *Q, *c, *a1, *P1, *P1_inf;
    double tol0, tolP;
    int inv_method;
    
    /*** Output variables ***/
    double *snlogL, *fvar;
    
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
    
    /*** Retrieve analysis settings ***/
    tol0        = mxGetScalar(prhs[15]);
    tolP        = mxGetScalar(prhs[16]);
    inv_method  = (int)mxGetScalar(prhs[17]);
    
    /*** Truncate Z and H w.r.t. missing data ***/
    if(miss) truncate_yHZ(p, m, n, N, &y, miss, &pmiss, &imiss, &H, &Hd, &Z, &Zd);
    else { pmiss = (int *)0; imiss = (int *)0; }

    /*** Allocate output variables ***/
    snlogL  = mxGetPr(plhs[0] = mxCreateDoubleMatrix(1, N, mxREAL));
    fvar    = (nlhs > 1) ? mxGetPr(plhs[1] = mxCreateDoubleMatrix(1, N, mxREAL)) : (double *)0;
    
    /*** Kalman filter ***/
    kalman(p, m, r, n, N, y, (int *)0, false, pmiss, sta, H, Hd, Z, Zd, T, Td, R, Rd, Q, Qd, c, cd, a1, P1, P1_inf,
        (double *)0, (double *)0, (double *)0, (int *)0, (int *)0, (double *)0, (double *)0, (double *)0, (double *)0,
        (double *)0, (double *)0, (double *)0, (double *)0, snlogL, fvar, (double *)0, (double *)0, (double *)0,
        tol0, tolP, inv_method);

    /*** Cleanup ***/
    cleanup();
}
