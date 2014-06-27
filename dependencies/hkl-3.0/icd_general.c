#include "mex.h"
#include <math.h>
    
 void mexFunction(int nlhs,
                              mxArray *plhs[],
                              int nrhs,
                              const mxArray *prhs[])
         
             {
double *z,a,b,c,maxdiagG;
double alpha,tol,*temp,*diagG,*G, *Gbis;
int m, n,i,j,jast;
int iter;
int *pp;
int nmax;
int nmax_to_allocate;
double *K, residual,trK,*diagK;

n = mxGetN(prhs[0]); /* number of samples */
K = mxGetPr(prhs[0]); 
temp=mxGetPr(prhs[1]);
tol=*temp;              /* approximation parameter */
temp=mxGetPr(prhs[2]);
	nmax=*temp;         /* maximal rank */


diagG= (double*) calloc (n,sizeof(double));
diagK= (double*) calloc (n,sizeof(double));
pp= (int*) calloc (n,sizeof(int));

plhs[0]=mxCreateDoubleMatrix(n,nmax,0); 
G= mxGetPr(plhs[0]); 


iter=0;

for (i=0;i<=n-1;i++)  pp[i]=i;
trK = 0;
for (i=0;i<=n-1;i++) 
{
     diagK[i]=K[i+n*i];
     diagG[i]=diagK[i];
     trK += diagK[i];
}
residual=trK;
jast=0;

while ( ( residual > tol * trK) & ( iter < nmax ) )
{

/* switches already calculated elements of G and order in pp */
if (jast!=iter)
	{
	i=pp[jast];  pp[jast]=pp[iter];  pp[iter]=i;
	for (i=0;i<=iter;i++)
		{
		a=G[jast+n*i];  G[jast+n*i]=G[iter+n*i];  G[iter+n*i]=a;
		}
	}


G[iter*(n+1)]=sqrt(diagG[jast]);
a=-alpha;

for (i=iter+1; i<=n-1; i++) 
	{
	G[i+n*iter]=K[pp[iter]+n*pp[i]];
    }

if (iter>0)
	for (j=0; j<=iter-1; j++)
		for (i=iter+1; i<=n-1; i++) G[i+n*iter]-=G[i+n*j]*G[iter+n*j];

for (i=iter+1; i<=n-1; i++) 
	{
	G[i+n*iter]/=G[iter*(n+1)];

	}
residual=0.0;
jast=iter+1;
maxdiagG=0;
for (i=iter+1; i<=n-1; i++)
	{
	b=diagK[pp[i]];
	for (j=0;j<=iter;j++)
		{
		 b-=G[i+j*n]*G[i+j*n];
		}
      diagG[i]=b;
	if (b>maxdiagG)
		{
		jast=i;
		maxdiagG=b;
		}
      residual+=b;
	} 

iter++;
}


plhs[1]=mxCreateDoubleMatrix(1,n,0); 
z= mxGetPr(plhs[1]); 
for (i=0;i<=n-1;i++) z[i]=1.0+pp[i];

plhs[2]=mxCreateDoubleMatrix(1,1,0); 
z= mxGetPr(plhs[2]); 
*z = iter;

plhs[3]=mxCreateDoubleMatrix(1,1,0); 
z= mxGetPr(plhs[3]); 
*z = residual / trK ;

free(diagG);
free(diagK);

free(pp);

             }


