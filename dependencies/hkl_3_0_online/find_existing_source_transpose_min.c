#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    int i, j, p, n, d, k;
    int *t, *T;
    double *aa, *bb, a, b;
    t = mxGetPr(prhs[0]);       /* row of integers  to look for */
    T = mxGetPr(prhs[1]);       /* rows of integers  to look i */
    p =  mxGetM(prhs[1]);
    n =  mxGetN(prhs[1]);
    k = 0;
    a = p;
    b = 0;
    for (i=0;i<n;i++) {
        d = 0;
        for (j=0;j<p;j++) {
            
            if (t[j]!=T[k++])
                d++;
            
            if (d >= a) { /* stop if more than before */
                k = k + p - 1 - j;
                j = p;
            }
        }
        if (d < a) {
            a = d;
            b = i+1;
        }
        
        if (a==0)   /* stop if exact match */
            i=n;
    }
    
    
    plhs[0]=mxCreateDoubleMatrix(1, 1, 0);
    aa = mxGetPr(plhs[0]);
    plhs[1]=mxCreateDoubleMatrix(1, 1, 0);
    bb = mxGetPr(plhs[1]);
    aa[0]=a;
    bb[0]=b;
    return;
    
}