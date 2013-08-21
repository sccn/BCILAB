#include "mex.h"
#include <math.h>
    
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
int k,k1,i,j,n,ni;
float *K, *z;

K = mxGetPr(prhs[0]); 
n =  mxGetM(prhs[0]);
n = floor( sqrt( 2 * n ) );


plhs[0]=mxCreateNumericMatrix(n,n,mxSINGLE_CLASS,0); 
z= mxGetPr(plhs[0]); 
k=0;
k1=0;
for (j=0;j<=n-1;j++)
	{
	ni=j;	
	for (i=0;i<=j;i++)
		{
		z[i+k1]=K[k];
		z[ni]=K[k];
		k++;
		ni+=n;		
		}
	k1 += n;
	}
}


