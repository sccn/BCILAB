#include "ssa_matlab.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*****************************************************************************\
     MATLAB
     [a P v F] = kalman_int(y, H, Hd, Z, Zd, T, Td, R, Rd, Q, Qd, c, cd, a1, P1,
                            hideinit, tol0, tolP, inv_method, postproc)
    \*****************************************************************************/
    
    /*** Input variables ***/
    int p, m, r, n, N;
    double *y;
    int *miss, *pmiss, *imiss;
    int Hd, Zd, Td, Rd, Qd, cd, sta;
    double *H, *Z, *T, *R, *Q, *c, *a1, *P1, *P1_inf;
    double tol0, tolP;
    int hideinit, inv_method, postproc;
    
    /*** Output variables ***/
    int dims[3];
    double *a, *P, *v, *F;
    int d;
    
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
    sta         = !Hd && !Zd && !Td && !Rd && !Qd;
    
    /*** Retrieve analysis settings ***/
    hideinit    = (int)mxGetScalar(prhs[15]);
    tol0        = mxGetScalar(prhs[16]);
    tolP        = mxGetScalar(prhs[17]);
    inv_method  = (int)mxGetScalar(prhs[18]);
    postproc    = (int)mxGetScalar(prhs[19]);
    
    /*** Truncate Z and H w.r.t. missing data ***/
    if(miss) truncate_yHZ(p, m, n, N, &y, miss, &pmiss, &imiss, &H, &Hd, &Z, &Zd);
    else { pmiss = (int *)0; imiss = (int *)0; }
    
    /*** Allocate output variables ***/
    dims[0] = m; dims[1] = N; dims[2] = n+1;
    a       = mxGetPr(plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
    dims[0] = m; dims[1] = m; dims[2] = n+1;
    P       = mxGetPr(plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
    if(miss) {
        v       = (double *)mxCalloc(p*N*n, sizeof(double));
        F       = (double *)mxCalloc(p*p*n, sizeof(double));
    }
    else {
        dims[0] = p; dims[1] = N; dims[2] = n;
        v       = mxGetPr(plhs[2] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
        dims[0] = p; dims[1] = p; dims[2] = n;
        F       = mxGetPr(plhs[3] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
    }
    
    /*** Kalman filter ***/
    kalman(p, m, r, n, N, y, (int *)0, false, pmiss, sta, H, Hd, Z, Zd, T, Td, R, Rd, Q, Qd, c, cd, a1, P1, P1_inf,
        a, P, (double *)0, &d, (int *)0, v, F, (double *)0, (double *)0, (double *)0, (double *)0, (double *)0, (double *)0,
        (double *)0, (double *)0, (double *)0, (double *)0, (double *)0, tol0, tolP, inv_method);
    
    /*** Postprocessing ***/
    if(miss) {
        double *vtemp, *Ftemp;
        dims[0] = p; dims[1] = N; dims[2] = n;
        vtemp   = mxGetPr(plhs[2] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
        dims[0] = p; dims[1] = p; dims[2] = n;
        Ftemp   = mxGetPr(plhs[3] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
        
        untruncate(p, n, N, (int *)0, false, pmiss, imiss, v, F, vtemp, Ftemp);
        
        mxFree(v);
        mxFree(F);
    }
    if(hideinit && d > 0) {
        int i;
        double NaN  = mxGetNaN();
        for(i=0; i<m*N*d; i++) a[i] = NaN;
        for(i=0; i<m*m*d; i++) P[i] = NaN;
        for(i=0; i<p*N*d; i++) v[i] = NaN;
        for(i=0; i<p*p*d; i++) F[i] = NaN;
    }
    if(postproc) {
        permute_data(plhs[0], m, N, n+1);
        permute_data(plhs[2], p, N, n);
        if(p == 1) permute_data(plhs[3], 1, 1, n);
        if(m == 1) permute_data(plhs[1], 1, 1, n+1);
    }
    
    /*** Cleanup ***/
    cleanup();
}
