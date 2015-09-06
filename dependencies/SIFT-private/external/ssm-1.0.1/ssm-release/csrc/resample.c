#include "mex.h"

void resample(double *y, int n, const double *w, const double *x)
{
    int i, j;
    double *cumw = (double *)mxCalloc(n, sizeof(double));

    cumw[0]     = w[0];
    for(i=1; i<n; i++)
        cumw[i] = cumw[i-1] + w[i];

    for(i=0; i<n; i++) {
        j       = 0;
        while(x[i] > cumw[j]) j++;
        y[i]    = j+1;
    }
    
    mxFree(cumw);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int m       = mxGetM(prhs[0]);
    int n       = mxGetN(prhs[0]);
    double *w   = mxGetPr(prhs[0]);
    double *x   = mxGetPr(prhs[1]);
    double *y;

    n           = (m > n) ? m : n;
    plhs[0]     = mxCreateDoubleMatrix(1, n, mxREAL);
    y           = mxGetPr(plhs[0]);
    
    resample(y, n, w, x);
}
