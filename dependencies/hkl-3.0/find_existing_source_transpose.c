#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i,j,p,n,d,k;
    int *t, *T;
    double *z;
    t = mxGetPr(prhs[0]);       /* row of integers  to look for */
    T = mxGetPr(prhs[1]);       /* rows of integers  to look i */
    p =  mxGetM(prhs[1]);
    n =  mxGetN(prhs[1]);
    plhs[0]=mxCreateDoubleMatrix(n,1,0);
    z= mxGetPr(plhs[0]);
    k = 0;
    for (i=0;i<n;i++)
    {
        d = 0;
        for (j=0;j<p;j++)
        { 
            /* mexPrintf("i=%d j=%d t=%d T=%d\n",i,j,t[j],T[i+n*j]); */
            if (t[j]!=T[k++])
                d++;
        }
        z[i] = d;
    }
    
}
