#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
int k,k1,i,j,n;
double  *z, *alpha,res;
float *K;
K = mxGetPr(prhs[0]);
alpha = mxGetPr(prhs[1]);
n =  mxGetM(prhs[1]);

res=0;
k=0;
k1=0;
for (j=0;j<=n-1;j++)
       {
       for (i=0;i<=j;i++)
               {
               if (i!=j) res+= 2 * ( ( double) K[k] )*alpha[i]*alpha[j] ;
               else
                   res+= ( ( double) K[k] )*alpha[j]*alpha[i];
               k++;
               }
       k1 += n;
       }


plhs[0]=mxCreateDoubleMatrix(1,1,0);
z= mxGetPr(plhs[0]);
z[0] = res;
}