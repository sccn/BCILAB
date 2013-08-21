#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i,j,p,n,nhull,k,s,ii,jj,iii,jjj,isok,delta,nsource,pi;
    int *hull, *qs,t, *children, nchildren;
    double  *output;
    hull = mxGetPr(prhs[0]);
    nhull =  mxGetN(prhs[0]);
    p =  mxGetM(prhs[0]);
    qs = mxGetPr(prhs[1]);
    
    
    
/* mexPrintf("p=%d - nhull=%d\n\n",p,nhull); */
        /* mexPrintf("p=%d - nhull=%d\n\n",p,nhull); */

    
    /* we may assume that all variables are active */
    /* consider all children which are not in the active set and for which all parents are in the active set */
    children = (int*) calloc (nhull*p*p,sizeof(int));
    k = 0;
/*     mexPrintf("p=%d - nhull=%d\n\n",p,nhull); */

    for (i=0;i<nhull;i++)                   /* loop over member of the hull */
        for (j=0;j<p;j++)             /* loop over active variables */
            if ( (  hull[p*i+j] ) < qs[j] )        /* check that it is not at the end of the line */
    {
        pi = p*i;
      /*  mexPrintf("%d \n",i); */
       /* for (s=0;s<p;s++)
                        if (s==j) mexPrintf("%d ", (  hull[p*i+indices[s]] )+1);
                        else  mexPrintf("%d ",  hull[p*i+indices[s]]);
       */
        
        
        /* first check that the children obtained by following the j-th direction is not already in the list of sources */
        isok = 1;
        ii=0;
        while ((isok>0) & (ii<k))
        {
            delta = 0;
            for (jj=0;jj<p;jj++)
            {
                if (jj==j) delta = abs( children[jj+ii*p] -  ( hull[pi+jj] ) - 1 );
                else delta = abs( children[jj+ii*p]-  ( hull[pi+jj] ));
                if (delta>0) jj = p;   /* gets out when sure */
            }
            if (delta==0) isok = 0;
            ii++;
        }
        
        
        
        /* check that the children obtained by following the j-th direction is not already in the hull */
        if (isok>0)
        {
            /* mexPrintf("a"); */
            ii=0;
            while ((isok>0) & (ii<nhull))
            {
                delta = 0;
                for (jj=0;jj<p;jj++)
                {
                    if (jj==j) delta = abs( (  hull[p*ii+jj] ) -  (  hull[pi+jj] ) - 1 );
                    else delta = abs( (hull[p*ii+jj] ) -  (  hull[pi+jj] ));
                    if (delta>0) jj = p;  /* gets out when sure */
                }
                if (delta==0) isok = 0;
                ii++;
            }
            
            if (isok>0)
            {
                /* mexPrintf("b"); */
                
            /* not in the active set -> check parents: they should all be in the active set */
                isok = 1;
                jj = 0;                         /* jj = index of parents of the candidate */
                
                while ((isok>0) & (jj<p))
                {
                    if ( (jj!=j) & hull[pi+jj] > 1 ) /* only considers the ones which are not at the origin */
                    {
                        isok=0;
                        for (ii=0;ii<nhull;ii++)
                        {
                        /* check equality with the members of the hull */
                            
                            delta = 0;
                            for (jjj=0;jjj<p;jjj++)       /* index on parent candidate */
                            {
                                if (jjj==jj)
                                    delta += abs( (  hull[p*ii+jj] ) -  (  hull[pi+jjj] ) + 1 );
                                else if (jjj==j)
                                    delta += abs( (  hull[p*ii+jjj] ) - 1 -  (  hull[pi+jjj] ));
                                else
                                    delta += abs( (  hull[p*ii+jjj] ) -  (  hull[pi+jjj] ));
                                if (delta>0) jjj=p; /* gets out when sure */
                            }
                            if (delta==0) { isok = 1; ii=nhull; }
                          /*  mexPrintf(" ii=%d jj=%d d=%d ",ii,jj,delta); */
                        }
                    }
                 /*   mexPrintf(" isok=%d ", isok); */
                    jj++;
                }
                if (isok>0)
                {
                    /* mexPrintf("c\n");*/
                    for (s=0;s<p;s++)
                        if (s==j) children[k*p+s] =  (  hull[pi+s] )+1;
                        else children[k*p+s] =  hull[pi+s];
                 /*   mexPrintf("*"); */
                    
                    k++;
                }
            }
        }
            }
    nchildren = k;
    nsource = nchildren;
    
    
    plhs[0]=mxCreateDoubleMatrix(nsource,p,0);
    output = mxGetPr(plhs[0]);
    for (i=0;i<nsource*p;i++) output[i]=1;
    
    
    k = 0;
    for (i=0;i<nchildren;i++)
    {
        for (s=0;s<p;s++)
            output[k+s*nsource ] = children[i*p+s];
        k++;
    }
    
    
    free(children);
    return;
}