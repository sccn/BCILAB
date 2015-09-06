#include "ssa_matlab.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /***************************************************************************************\
     MATLAB
     [alphahat V r N] = statesmo_int(y, H, Hd, Z, Zd, T, Td, R, Rd, Q, Qd, c, cd, a1, P1,
                                        tol0, tolP, inv_method, postproc)
    \***************************************************************************************/
    
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
    double *alphahat, *V, *r, *N;
    
    /*** Buffer variables ***/
    double *a, *P, *P_inf, *v, *invF, *F2, *L, *L1;
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
    dims[0]     = m; dims[1] = NN; dims[2] = n;
    alphahat    = mxGetPr(plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
    dims[0]     = m; dims[1] = m; dims[2] = n;
    V           = mxGetPr(plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
    dims[0]     = m; dims[1] = NN; dims[2] = n+1;
    r           = mxGetPr(plhs[2] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
    dims[0]     = m; dims[1] = m; dims[2] = n+1;
    N           = mxGetPr(plhs[3] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
    
    /*** Allocate workspace ***/
    a       = (double *)mxCalloc(m*NN*(n+1), sizeof(double));
    P       = (double *)mxCalloc(m*m*(n+1), sizeof(double));
    P_inf   = (double *)mxCalloc(m*m*(n+1), sizeof(double));
    Fns     = (int *)mxCalloc(n, sizeof(int));
    v       = (double *)mxCalloc(p*NN*n, sizeof(double));
    invF    = (double *)mxCalloc(p*p*n, sizeof(double));
    F2      = (double *)mxCalloc(p*p*n, sizeof(double));
    L       = (double *)mxCalloc(m*m*n, sizeof(double));
    L1      = (double *)mxCalloc(m*m*n, sizeof(double));

    /*** Kalman filter ***/
    kalman(p, m, rr, n, NN, y, (int *)0, false, pmiss, sta, H, Hd, Z, Zd, T, Td, R, Rd, Q, Qd, c, cd, a1, P1, P1_inf,
        a, P, P_inf, &d, Fns, v, (double *)0, invF, F2, (double *)0, (double *)0,
        L, L1, (double *)0, (double *)0, (double *)0, (double *)0, (double *)0, tol0, tolP, inv_method);
    
    /*** Fast smoother ***/
    statesmo(p, m, n, NN, (int *)0, false, pmiss, Z, Zd, T, Td, a, P, P_inf, d, Fns, v, invF, F2, L, L1,
        r, (double *)0, N, (double *)0, (double *)0, alphahat, V);
    
    /*** Cleanup ***/
    cleanup();
    mxFree(a);
    mxFree(P);
    mxFree(P_inf);
    mxFree(Fns);
    mxFree(v);
    mxFree(invF);
    mxFree(F2);
    mxFree(L);
    mxFree(L1);
    
    /*** Postprocessing ***/
    if(postproc) {
        permute_data(plhs[0], m, NN, n);
        permute_data(plhs[2], m, NN, n+1);
        if(m == 1) {
            permute_data(plhs[1], 1, 1, n);
            permute_data(plhs[3], 1, 1, n+1);
        }
    }
}
