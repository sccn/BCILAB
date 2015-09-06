#include "mex.h"

void diaginprod(double *y, int m, int n, const double *A, const double *x)
{
    int i, j, k;
    
    for(i=0; i<n; i++) {
        int im  = i*m;
        y[i]    = 0;
        for(j=0; j<m; j++) {
            double Ax   = 0;
            for(k=0; k<m; k++)
                Ax      = Ax + A[j+k*m]*x[k+im];
            y[i]        = y[i] + x[j+im]*Ax;
        }
    }
}

void diffdiaginprod(double *y, int m, int n, const double *A, const double *x1, const double *x2)
{
    /* A is a m x m square matrix
       x1, x2 is m x n matrices representing n column vectors
       y[i, j] = (x1[j]-x2[i])'*A*(x1[j]-x2[i]) for each column x1[j], x2[i]
    */
    int i, j, k, l;
    double *x = (double *)mxCalloc(m, sizeof(double));
    
    for(i=0; i<n; i++) {
        for(j=0; j<n; j++) {
            int jm      = j*m;
            int ijn     = i + j*n;
            
            for(k=0; k<m; k++)
                x[k]    = x1[k+jm] - x2[k+i*m];
            y[ijn]      = 0;
            for(k=0; k<m; k++) {
                double Ax   = 0;
                for(l=0; l<m; l++)
                    Ax      = Ax + A[k+l*m]*x[l];
                y[ijn]      = y[ijn] + x[k]*Ax;
            }
        }
    }
    
    mxFree(x);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int m       = mxGetM(prhs[0]);
    int n       = mxGetN(prhs[1]);
    double *A   = mxGetPr(prhs[0]);
    double *y;
    
    if(nrhs == 2) {
        double *x   = mxGetPr(prhs[1]);
        
        plhs[0]     = mxCreateDoubleMatrix(n, 1, mxREAL);
        y           = mxGetPr(plhs[0]);
        
        diaginprod(y, m, n, A, x);
    }
    else {
        double *x1  = mxGetPr(prhs[1]);
        double *x2  = mxGetPr(prhs[2]);
        
        plhs[0]     = mxCreateDoubleMatrix(n, n, mxREAL);
        y           = mxGetPr(plhs[0]);
        
        diffdiaginprod(y, m, n, A, x1, x2);
    }
}
