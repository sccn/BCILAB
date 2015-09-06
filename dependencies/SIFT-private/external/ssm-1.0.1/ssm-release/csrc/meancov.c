#include "mex.h"

void meancov(double *mean, double *cov, int m, int n, int N, const double *x)
{
    int t, mn   = m*n;
    double *x2  = (double *)mxCalloc(m*N, sizeof(double));
    
    for(t=0; t<n; t++) {
        int tm  = t*m;
        int tmm = tm*m;
        int i, j, k;

        for(i=0; i<m; i++) {
            int itm = i + tm;
            mean[itm]   = 0;
            for(j=0; j<N; j++)
                mean[itm]   += x[itm+j*mn];
            mean[itm]   /= N;
            for(j=0; j<N; j++)
                x2[i+j*m]   = x[itm+j*mn] - mean[itm];
        }
        
        for(i=0; i<m; i++) {
            int itmm    = i + tmm;
            for(j=i; j<m; j++) {
                int ijtmm       = itmm + j*m;
                cov[ijtmm]      = 0;
                for(k=0; k<N; k++)
                    cov[ijtmm]  += x2[i+k*m]*x2[j+k*m];
                cov[ijtmm]      /= N-1;
                cov[j+i*m+tmm]  = cov[ijtmm];
            }
        }
    }

    mxFree(x2);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const int *dim  = mxGetDimensions(prhs[0]);
    int m           = dim[0];
    int n           = dim[1];
    int N           = dim[2];
    int covdim[3];
    double *x       = mxGetPr(prhs[0]);
    double *mean;
    double *cov;

    covdim[0]   = m;
    covdim[1]   = m;
    covdim[2]   = n;
    plhs[0]     = mxCreateDoubleMatrix(m, n, mxREAL);
    plhs[1]     = mxCreateNumericArray(3, covdim, mxDOUBLE_CLASS, mxREAL);
    mean        = mxGetPr(plhs[0]);
    cov         = mxGetPr(plhs[1]);
    
    meancov(mean, cov, m, n, N, x);
}
