#include "ssa_matlab.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /**************************************************************************************************\
     MATLAB
     [y alpha eps eta] = sample_int(n, N, H, Hd, Z, Zd, T, Td, R, Rd, Q, Qd, c, cd, a1, P1, postproc)
     sample_int('seed', S)
    \**************************************************************************************************/
    
    if(mxIsChar(prhs[0])) {
        int length          = mxGetNumberOfElements(prhs[0]) + 1;
        char *str           = (char *)mxCalloc(length, sizeof(char));
        
        mxGetString(prhs[0], str, length);
        if(strcmp(str, "seed") == 0) {
            int n               = mxGetNumberOfElements(prhs[1]);
            unsigned long *seed = (unsigned long *)mxGetData(prhs[1]);
            if(n == 1) init_genrand(*seed);
            else init_by_array(seed, n);
        }
        
        mxFree(str);
        return;
    }
    else {
        /*** Input variables ***/
        int p, m, r, n, N;
        int Hd, Zd, Td, Rd, Qd, cd;
        double *H, *Z, *T, *R, *Q, *c, *a1, *P1, *P1_inf;
        int postproc;
        
        /*** Output variables ***/
        int dims[3];
        double *y, *alpha, *eps, *eta;
        
        /*** Buffer variables ***/
        double *samples;
        
        /*** Retrieve data ***/
        n   = (int)mxGetScalar(prhs[0]);
        N   = (int)mxGetScalar(prhs[1]);
        
        /*** Retrieve model ***/
        get_M(prhs + 2, &p, &H, &Hd);
        get_M(prhs + 4, (int *)0, &Z, &Zd);
        get_M(prhs + 6, &m, &T, &Td);
        get_M(prhs + 8, (int *)0, &R, &Rd);
        get_M(prhs + 10, &r, &Q, &Qd);
        get_M(prhs + 12, (int *)0, &c, &cd);
        get_M(prhs + 14, (int *)0, &a1, (int *)0);
        get_P1(prhs[15], (int *)0, &P1, &P1_inf);
        
        /*** Retrieve analysis settings ***/
        postproc    = (int)mxGetScalar(prhs[16]);
        
        /*** Generate samples ***/
        samples = (double *)mxCalloc((m+(p+r)*n)*N, sizeof(double));
        genrandn((m+(p+r)*n)*N, samples);
        
        /*** Allocate output variables ***/
        dims[0] = p; dims[1] = N; dims[2] = n;
        y       = mxGetPr(plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
        dims[0] = m; dims[1] = N; dims[2] = n;
        alpha   = mxGetPr(plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
        dims[0] = p; dims[1] = N; dims[2] = n;
        eps     = mxGetPr(plhs[2] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
        dims[0] = r; dims[1] = N; dims[2] = n;
        eta     = mxGetPr(plhs[3] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL));
        
        /*** Sample ***/
        sample(p, m, r, n, N, H, Hd, Z, Zd, T, Td, R, Rd, Q, Qd, c, cd, a1, P1, samples, y, alpha, eps, eta);
        
        /*** Cleanup ***/
        mxFree(samples);
        cleanup();
        
        /*** Postprocessing ***/
        if(postproc) {
            permute_data(plhs[0], p, N, n);
            permute_data(plhs[1], m, N, n);
            permute_data(plhs[2], p, N, n);
            permute_data(plhs[3], r, N, n);
        }
    }
}
