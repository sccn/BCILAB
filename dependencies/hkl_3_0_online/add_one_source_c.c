#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    int *rest_of_sources, nrest, *qs, *toadd, *added;
    int k, j, i, p, s, pi, jj, isok, isdes;
    double *output;
    rest_of_sources = mxGetPr(prhs[0]);
    nrest =  mxGetN(prhs[0]);
    p =  mxGetM(prhs[0]);
    toadd = mxGetPr(prhs[1]);
    qs = mxGetPr(prhs[2]);
    
    
    /*mexPrintf("p=%d \n",p); */
    
    /* consider all children of source to add */
    added = (int*) calloc(p*p, sizeof(int));
    k = 0;
    
    
    for (j=0;j<p;j++)             /* loop over children */
        if ( (  toadd[j] ) < qs[j] )        /* check that it is not at the end of the line */ {
        
        /* loop over rest of sources, as soon as it is the descendent of one of them, stop and go to the next */
        isok = 1;
        i=0;
        while ((isok>0) & (i<nrest)) {
            pi = p*i;
            isdes = 1;
            for (jj=0;jj<p;jj++) {
                if (jj==j) {
                    if (rest_of_sources[jj+pi] >  toadd[jj] + 1 ) {
                        isdes = 0;      /* the particular other source is not a descendant */
                        jj=p;
                    }
                }
                else {
                    if (rest_of_sources[jj+pi] >  toadd[jj] ) {
                        isdes = 0;      /* the particular other source is not a descendant */
                        jj=p;
                    }
                }
            }
            if (isdes>0)
                isok = 0;      /* stop because we have a descendant*/
            i++;
        }
        
        if (isok>0) {
            for (s=0;s<p;s++)
                if (s==j) added[k*p+s] = toadd[s] +1;
                else added[k*p+s] =   toadd[s];
            k++;
        }
        
        }
    
    
    plhs[0]=mxCreateDoubleMatrix(k, p, 0);
    output = mxGetPr(plhs[0]);
    
    for (i=0;i<k;i++)
        for (s=0;s<p;s++)
            output[i+s*k ] = added[i*p+s];
    
    
    free(added);
    return;
}