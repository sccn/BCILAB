#include "mex.h"
#include <math.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i,j,p,n,nhull,k,s,ii,jj,iii,jjj,isok,delta,nsource,pi;
    int *hull, *qs,t, *active,*indices,nactive, *children, nchildren;
    double  *output;
    hull = mxGetPr(prhs[0]);
    nhull =  mxGetN(prhs[0]);
    p =  mxGetM(prhs[0]);
    qs = mxGetPr(prhs[1]);
    
    
    /* first find out which variables are used, which are not */
    active = (int*) calloc (p,sizeof(int));
    indices = (int*) calloc (p,sizeof(int));
    for (i=0;i<p;i++) active[i] =0;
    k=0;
    for (i=0;i<nhull;i++)
           for (j=0;j<p;j++)
     if  (hull[k++]>1) active[ j ] = 1;
    k=0;
    for (j=0;j<p;j++) if (active[j]==1) indices[k++] = j;
    nactive = k;
    
/*
 for (i=0;i<p;i++) mexPrintf("%d ",active[i]);
    mexPrintf("\n");
    for (i=0;i<p;i++) mexPrintf("%d ",indices[i]);
 */
    
    
    /* now consider only the ones which are active */
    /* consider all children which are not in the active set and for which all parents are in the active set */
    children = (int*) calloc (nhull*nactive*nactive,sizeof(int));
    k = 0;
    
    for (i=0;i<nhull;i++)                   /* loop over member of the hull */
        for (j=0;j<nactive;j++)             /* loop over active variables */
            if ( (  hull[p*i+indices[j]] ) < qs[j] )        /* check that it is not at the end of the line */
    {
      /*  mexPrintf("\n");
		for (s=0;s<nactive;s++)
                        if (s==j) mexPrintf("%d ", (  hull[p*i+indices[s]] )+1);
                        else  mexPrintf("%d ",  hull[p*i+indices[s]]);
		*/
		
      
        /* first check that the children obtained by following the j-th direction is not already in the list of sources */
        isok = 1;
        ii=0;
        while ((isok>0) & (ii<k))
        {
            delta = 0;
            for (jj=0;jj<nactive;jj++)
            {
                if (jj==j) delta = abs( children[jj+ii*nactive] -  ( hull[p*i+indices[jj]] ) - 1 );
                else delta = abs( children[jj+ii*nactive]-  ( hull[p*i+indices[jj]] ));
                if (delta>0) jj = nactive;   /* gets out when sure */
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
                for (jj=0;jj<nactive;jj++)
                {
                    if (jj==j) delta = abs( (  hull[p*ii+indices[jj]] ) -  (  hull[p*i+indices[jj]] ) - 1 );
                    else delta = abs( (hull[p*ii+indices[jj]] ) -  (  hull[p*i+indices[jj]] ));
                    if (delta>0) jj = nactive;  /* gets out when sure */
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
                
                while ((isok>0) & (jj<nactive))
                {
                    if ( (jj!=j) & hull[p*i+indices[jj]] > 1 ) /* only considers the ones which are not at the origin */
                    {
                        isok=0;
                        for (ii=0;ii<nhull;ii++)
                        {
                        /* check equality with the members of the hull */
                            
                            delta = 0;
                            for (jjj=0;jjj<nactive;jjj++)       /* index on parent candidate */
                            {
                                if (jjj==jj)
                                    delta += abs( (  hull[p*ii+indices[jjj]] ) -  (  hull[p*i+indices[jjj]] ) + 1 );
                                else if (jjj==j)
                                    delta += abs( (  hull[p*ii+indices[jjj]] ) - 1 -  (  hull[p*i+indices[jjj]] ));
                                else
                                    delta += abs( (  hull[p*ii+indices[jjj]] ) -  (  hull[p*i+indices[jjj]] ));
                                if (delta>0) jjj=nactive; /* gets out when sure */
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
                    for (s=0;s<nactive;s++)
                        if (s==j) children[k*nactive+s] =  (  hull[p*i+indices[s]] )+1;
                        else children[k*nactive+s] =  hull[p*i+indices[s]];
                 /*   mexPrintf("*"); */
                    
                    k++;
                }
            }
        }
            }
    nchildren = k;
    nsource = nchildren+p-nactive;
    
    
    plhs[0]=mxCreateDoubleMatrix(nsource,p,0);
    output = mxGetPr(plhs[0]);
    for (i=0;i<nsource*p;i++) output[i]=1;
    
    
    k = 0;
    for (i=0;i<nchildren;i++)
    {
        for (s=0;s<nactive;s++)
        {
            output[k+indices[s]*nsource ] = children[i*nactive+s];
        /*     mexPrintf("*  %d %d %d ",k,indices[s],children[i*nactive+s]); */
            
        }
        k++;
    }
    
    for (j=0;j<p;j++)
        if (active[j]==0)
            output[k++ + j*nsource ] = 2;
    
    
    
    free(active);
    free(indices);
    free(children);
    return;
}
