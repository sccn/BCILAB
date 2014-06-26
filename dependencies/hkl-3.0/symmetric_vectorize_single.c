#include "mex.h"
#include <math.h>
    
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
int k,k1,i,j,n;
float *K, *z;

K = mxGetPr(prhs[0]); 
n =  mxGetM(prhs[0]);

plhs[0]=mxCreateNumericMatrix(n*(n+1)/2,1,mxSINGLE_CLASS,0); 
z= mxGetPr(plhs[0]); 
k=0;
k1=0;
for (j=0;j<=n-1;j++)
	{	
	for (i=0;i<=j;i++)
		{
		z[k]=K[i+k1];
		k++;		
		}
	k1 += n;
	}
}


